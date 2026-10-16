"""Aim 1, step 3: a deployable wrist sleep model - nothing may depend on the Dreem recording.

Part 1  continuous 30-s feature table per participant (phase0/features30s/<ID>.parquet)
  grid: every 30 s from first to last wrist sample (UTC); local time from the participant tz.
  'night' = local noon-to-noon window, labelled by the date of the evening.
  worn: accelerometer data present and the 60-min SD rule of 04_daily_qc not triggered
        (recomputed here on 30-s means for speed) and vendor wear % >= 50 where available.
  features: ENMO mean/max, log ENMO, max 5-s angle change, mean anglez, HR, RMSSD, n beats,
            temperature, EDA; HR/RMSSD/temp/EDA z-scored WITHIN the noon-noon night (worn epochs);
            clock time (sin/cos of local hour); centred rolling means (1, 2.5, 5, 10, 30, 60 min)
            and still fraction (angle change <= 5 deg).  No vendor algorithm outputs.
Part 2  model: XGBoost binary sleep/wake and 4-class (W/light/deep/REM) trained on Dreem-labelled
  epochs (QC-pass recordings, >= 80 % wrist coverage, MT/U excluded), 5 folds grouped by family.
  Out-of-fold probabilities are kept for labelled epochs; a final model fitted on all labelled
  epochs scores every other epoch (unlabelled nights only - never a labelled one).
Part 3  sleep period time (SPT) per night, two methods:
  model : rolling 30-min median of P(sleep) >= 0.5 on worn epochs; runs separated by < 60 min
          merged; longest run = SPT (min 60 min).
  HDCZA : van Hees et al. 2018 - 5-s |d anglez|, 5-min rolling median, threshold = 15 x P10 of
          that series within the night, candidate blocks > 30 min, gaps < 60 min merged, longest = SPT.
  Evaluated on Dreem nights against Dreem sleep onset (first sleep epoch) and final awakening
  (last sleep epoch); TST = model-sleep epochs inside the model SPT vs Dreem TST.
Usage: python -I 12_aim1_deployable.py <phase0_dir> <aim1_dir> [n_workers]
"""
import os, sys, glob, json
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

P0, A1 = sys.argv[1:3]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
FT = os.path.join(P0, 'features30s'); os.makedirs(FT, exist_ok=True)
R = os.path.join(A1, 'deployable'); os.makedirs(R, exist_ok=True)
WINS = (2, 5, 10, 20, 60, 120)


def runs(b):
    """start/end (exclusive) index pairs of True runs."""
    d = np.diff(np.r_[0, b.astype(int), 0]); return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def longest_block(b, min_len, max_gap):
    s, e = runs(b)
    if not len(s):
        return None
    blocks = [[s[0], e[0]]]
    for a, z in zip(s[1:], e[1:]):
        if a - blocks[-1][1] < max_gap:
            blocks[-1][1] = z
        else:
            blocks.append([a, z])
    blocks = [bk for bk in blocks if bk[1] - bk[0] >= min_len]
    return max(blocks, key=lambda bk: bk[1] - bk[0]) if blocks else None


def build(pid):
    f5 = os.path.join(P0, 'epoch5s', f'{pid}.parquet')
    if not os.path.exists(f5):
        return pid, None
    if os.environ.get('REUSE_FEATURES') and os.path.exists(os.path.join(FT, f'{pid}.parquet')):
        return pid, 'reused'
    e5 = pd.read_parquet(f5)
    fm = os.path.join(P0, 'minute', f'{pid}.parquet')
    tz = pd.read_parquet(fm, columns=['tz']).tz.iloc[0] if os.path.exists(fm) else 'UTC'
    e5['dang'] = e5.anglez.diff().abs()
    e5['t30'] = e5.ts_utc.dt.floor('30s')
    g = e5.groupby('t30')
    f = pd.DataFrame({'enmo_mean': g.enmo.mean(), 'enmo_max': g.enmo.max(), 'dangle_max': g.dang.max(),
                      'anglez_mean': g.anglez.mean(), 'x': g.x.mean(), 'y': g.y.mean(), 'z': g.z.mean(),
                      'v_within': (g.sdx.mean() ** 2 + g.sdy.mean() ** 2 + g.sdz.mean() ** 2) / 3})
    grid = pd.date_range(f.index.min(), f.index.max(), freq='30s', tz='UTC')
    f = f.reindex(grid); f.index.name = 'ts_utc'
    f30 = os.path.join(P0, 'epoch30s', f'{pid}.parquet')
    if os.path.exists(f30):
        p = pd.read_parquet(f30).set_index('ts_utc')
        f = f.join(p, how='left')
    if os.path.exists(fm):
        m = pd.read_parquet(fm, columns=['ts_utc', 'wearing_detection_percentage']).set_index('ts_utc')
        f['vendor_wear'] = m.wearing_detection_percentage.reindex(f.index.floor('min')).values
    else:
        f['vendor_wear'] = np.nan
    # wear: 60-min windows (120 epochs) whose SD of 30-s axis means and within-epoch SD are both tiny
    has = f.enmo_mean.notna()
    still = pd.Series(False, index=f.index)
    for ax in 'xyz':
        sd = f[ax].rolling(120, center=True, min_periods=60).std()
        still = still.add((sd < 0.013) & (f.v_within.rolling(120, center=True, min_periods=60).mean() < 0.013 ** 2), fill_value=0)
    f['worn'] = has & (still < 2) & ~(f.vendor_wear < 50)
    loc = f.index.tz_convert(tz)
    f['ts_local'] = loc.tz_localize(None)
    hr = loc.hour + loc.minute / 60
    f['night'] = (loc - pd.Timedelta(hours=12)).date
    f['hour_sin'], f['hour_cos'] = np.sin(2 * np.pi * hr / 24), np.cos(2 * np.pi * hr / 24)
    f['log_enmo'] = np.log10(f.enmo_mean.clip(lower=1e-4))
    W = f.worn
    for c in ['enmo_mean', 'enmo_max', 'dangle_max', 'anglez_mean', 'log_enmo', 'hr_bpm', 'rmssd_ms', 'n_beats', 'temp_c', 'eda_us']:
        if c in f:
            f.loc[~W, c] = np.nan
    gn = f.groupby('night')
    for c in ['hr_bpm', 'rmssd_ms', 'temp_c', 'eda_us']:
        if c in f:
            f[c + '_z'] = (f[c] - gn[c].transform('median')) / gn[c].transform('std')
        else:
            f[c + '_z'] = np.nan
    f['still'] = (f.dangle_max <= 5).astype(float).where(f.dangle_max.notna())
    for w in WINS:
        for c in ['log_enmo', 'dangle_max', 'hr_bpm_z', 'rmssd_ms_z', 'temp_c_z', 'eda_us_z', 'still']:
            f[f'{c}_m{w}'] = f[c].rolling(w, center=True, min_periods=max(1, w // 3)).mean()
    f = f.drop(columns=['x', 'y', 'z', 'v_within'])
    f.reset_index().to_parquet(os.path.join(FT, f'{pid}.parquet'), index=False)
    return pid, dict(tz=tz, epochs=len(f), worn_frac=float(W.mean()))


def hdcza(pid):
    """van Hees 2018 SPT per noon-noon night from 5-s anglez."""
    e5 = pd.read_parquet(os.path.join(P0, 'epoch5s', f'{pid}.parquet'), columns=['ts_utc', 'anglez'])
    fm = os.path.join(P0, 'minute', f'{pid}.parquet')
    tz = pd.read_parquet(fm, columns=['tz']).tz.iloc[0] if os.path.exists(fm) else 'UTC'
    s = e5.set_index('ts_utc').anglez.resample('5s').mean()
    x = s.diff().abs().rolling(60, center=True, min_periods=30).median()        # 5 min
    loc = x.index.tz_convert(tz)
    night = pd.Series((loc - pd.Timedelta(hours=12)).date, index=x.index)
    out = []
    for n, xi in x.groupby(night):
        v = xi.values
        if np.isfinite(v).sum() < 12 * 60 * 6:                                  # need >= 6 h of data
            continue
        thr = 15 * np.nanpercentile(v, 10)
        b = longest_block(np.nan_to_num(v, nan=np.inf) < thr, 12 * 30, 12 * 60)
        if b:
            out.append(dict(empatica_id=pid, night=n, hdcza_on=xi.index[b[0]], hdcza_off=xi.index[min(b[1], len(v) - 1)]))
    return out


if __name__ == '__main__':
    L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
    ids = sorted(os.path.basename(f)[:-8] for f in glob.glob(os.path.join(P0, 'epoch5s', '*.parquet')))
    with ProcessPoolExecutor(NW) as ex:
        info = dict(ex.map(build, ids))
        hz = [r for rr in ex.map(hdcza, ids) for r in rr]
    print('feature build:', sum(v == 'reused' for v in info.values()), 'reused of', len(info))
    json.dump(info, open(os.path.join(FT, '_build_info.json'), 'w'), indent=1, default=str)
    HZ = pd.DataFrame(hz); HZ.to_csv(os.path.join(R, 'spt_hdcza.csv'), index=False)

    # ---- labelled epochs: Dreem epoch -> nearest 30-s grid point (<= 15 s) -------------------------
    import xgboost as xgb
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import cohen_kappa_score, f1_score
    D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
    pairs = pd.read_parquet(os.path.join(A1, 'aim1_epoch_pairs.parquet'), columns=['recording_id', 'empatica_id', 'coverage']).drop_duplicates('recording_id')
    keep = pairs[pairs.coverage >= 0.8]
    D = D.merge(keep[['recording_id', 'empatica_id']], on='recording_id')
    D = D[(D.qc_pass_core == True) & ~D.stage.isin(['MT', 'U'])].copy()
    D['ts_grid'] = D.ts_utc.dt.round('30s')
    lab = []
    for pid, d in D.groupby('empatica_id'):
        f = pd.read_parquet(os.path.join(FT, f'{pid}.parquet'))
        lab.append(d.merge(f, left_on='ts_grid', right_on='ts_utc', how='inner', suffixes=('_dreem', '')))
    T = pd.concat(lab, ignore_index=True)
    T = T.merge(L[['empatica_id', 'base_id', 'rarecnv_family_id', 'genotype']], on='empatica_id', how='left')
    T['family'] = T.rarecnv_family_id.fillna(T.base_id)
    T = T[T.worn]
    ndup = int(T.duplicated(['empatica_id', 'ts_grid']).sum())     # overlapping Dreem recordings
    T = T.drop_duplicates(['empatica_id', 'ts_grid'], keep='first').reset_index(drop=True)
    print('duplicate labelled epochs from overlapping Dreem recordings dropped:', ndup)
    T['sleep'] = T.stage.isin(['N1', 'N2', 'N3', 'R']).astype(int)
    T['y4'] = T.stage.map({'W': 0, 'N1': 1, 'N2': 1, 'N3': 2, 'R': 3})
    feats = [c for c in pd.read_parquet(os.path.join(FT, f'{ids[0]}.parquet')).columns
             if c not in ('ts_utc', 'ts_local', 'night', 'worn', 'vendor_wear', 'still')]
    X = T[feats].astype('float32').values
    prm = dict(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, tree_method='hist', n_jobs=NW)
    T['p_sleep'] = np.nan; T['pred4'] = np.nan
    for tr, te in GroupKFold(5).split(X, T.sleep, groups=T.family):
        T.loc[te, 'p_sleep'] = xgb.XGBClassifier(**prm, eval_metric='logloss').fit(X[tr], T.sleep.values[tr]).predict_proba(X[te])[:, 1]
        T.loc[te, 'pred4'] = xgb.XGBClassifier(**prm, objective='multi:softprob', num_class=4).fit(X[tr], T.y4.values[tr]).predict(X[te])
    yb = (T.p_sleep >= 0.5).astype(int)
    tn = ((yb == 0) & (T.sleep == 0)).sum(); fp = ((yb == 1) & (T.sleep == 0)).sum(); fn = ((yb == 0) & (T.sleep == 1)).sum(); tp = ((yb == 1) & (T.sleep == 1)).sum()
    f1 = f1_score(T.y4, T.pred4.astype(int), average=None, labels=[0, 1, 2, 3])
    EP = dict(epochs=len(T), recordings=int(T.recording_id.nunique()), participants=int(T.empatica_id.nunique()), families=int(T.family.nunique()),
              accuracy=float((tp + tn) / len(T)), sensitivity=float(tp / (tp + fn)), specificity=float(tn / (tn + fp)),
              kappa=float(cohen_kappa_score(T.sleep, yb)), acc4=float((T.y4 == T.pred4).mean()),
              kappa4=float(cohen_kappa_score(T.y4, T.pred4.astype(int))), F1_W_light_deep_REM=[round(float(v), 3) for v in f1])
    final_b = xgb.XGBClassifier(**prm, eval_metric='logloss').fit(X, T.sleep.values)
    final_4 = xgb.XGBClassifier(**prm, objective='multi:softprob', num_class=4).fit(X, T.y4.values)
    final_b.save_model(os.path.join(R, 'model_binary.json')); final_4.save_model(os.path.join(R, 'model_4class.json'))
    json.dump(feats, open(os.path.join(R, 'features.json'), 'w'))

    # ---- score every epoch, then SPT per night -----------------------------------------------------
    oof = T.set_index(['empatica_id', 'ts_grid'])[['p_sleep', 'pred4']]
    SP, NT = [], []
    for pid in ids:
        f = pd.read_parquet(os.path.join(FT, f'{pid}.parquet'))
        Xf = f[feats].astype('float32').values
        f['p_sleep'] = final_b.predict_proba(Xf)[:, 1]
        f['pred4'] = final_4.predict(Xf)
        if pid in oof.index.get_level_values(0):            # labelled epochs: use out-of-fold values
            o = oof.loc[pid]; ix = f.ts_utc.isin(o.index)
            f.loc[ix, ['p_sleep', 'pred4']] = o.reindex(f.loc[ix, 'ts_utc']).values
        f.loc[~f.worn, ['p_sleep', 'pred4']] = np.nan
        f[['ts_utc', 'night', 'worn', 'p_sleep', 'pred4']].to_parquet(os.path.join(R, f'scored_{pid}.parquet'), index=False)
        sm = f.p_sleep.rolling(61, center=True, min_periods=20).median()
        for n, idx in f.groupby('night').groups.items():
            v = sm.loc[idx].values
            if np.isfinite(v).sum() < 2 * 60 * 6:
                continue
            b = longest_block(np.nan_to_num(v, nan=0) >= 0.5, 120, 120)
            if b:
                seg = f.loc[idx[b[0]:b[1]]]
                SP.append(dict(empatica_id=pid, night=n, spt_on=seg.ts_utc.iloc[0], spt_off=seg.ts_utc.iloc[-1],
                               spt_min=len(seg) / 2, tst_min=float((seg.p_sleep >= 0.5).sum() / 2),
                               worn_frac=float(seg.worn.mean()),
                               **{f'min_{k}': float((seg.pred4 == i).sum() / 2) for i, k in enumerate(['W', 'light', 'deep', 'REM'])}))
    SP = pd.DataFrame(SP); SP.to_csv(os.path.join(R, 'spt_model_nightly.csv'), index=False)

    # ---- validate SPT on Dreem nights ---------------------------------------------------------------
    D2 = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
    D2 = D2.merge(keep[['recording_id', 'empatica_id']], on='recording_id')
    D2 = D2[D2.qc_pass_core == True]
    ref = D2.groupby('recording_id').apply(lambda d: pd.Series(dict(
        empatica_id=d.empatica_id.iloc[0], rec_start=d.ts_utc.min(), rec_end=d.ts_utc.max(),
        onset=d.ts_utc[d.stage.isin(['N1', 'N2', 'N3', 'R'])].min(), final_wake=d.ts_utc[d.stage.isin(['N1', 'N2', 'N3', 'R'])].max(),
        tst=d.stage.isin(['N1', 'N2', 'N3', 'R']).sum() / 2, n3=(d.stage == 'N3').sum() / 2, rem=(d.stage == 'R').sum() / 2,
        ends_in_wake=bool(d.stage[d.stage != 'MT'].iloc[-1] == 'W') if (d.stage != 'MT').any() else False)), include_groups=False).reset_index()
    out = []
    for r in ref.itertuples():
        for meth, tab, on, off in (('model', SP, 'spt_on', 'spt_off'), ('hdcza', HZ, 'hdcza_on', 'hdcza_off')):
            c = tab[tab.empatica_id == r.empatica_id]
            c = c[(c[on] < r.final_wake) & (c[off] > r.onset)]          # SPT overlapping the Dreem sleep
            if not len(c):
                out.append(dict(recording_id=r.recording_id, method=meth, found=False)); continue
            c = c.iloc[(np.minimum(c[off], r.final_wake) - np.maximum(c[on], r.onset)).values.argmax()]
            rec = dict(recording_id=r.recording_id, empatica_id=r.empatica_id, method=meth, found=True, ends_in_wake=r.ends_in_wake,
                       onset_err_min=(c[on] - r.onset).total_seconds() / 60, wake_err_min=(c[off] - r.final_wake).total_seconds() / 60,
                       ref_tst=r.tst, ref_n3=r.n3, ref_rem=r.rem)
            if meth == 'model':
                rec.update(est_tst=c.tst_min, est_deep=c.min_deep, est_rem=c.min_REM)
            out.append(rec)
    V = pd.DataFrame(out).merge(L[['empatica_id', 'genotype']], on='empatica_id', how='left')
    V.to_csv(os.path.join(R, 'spt_validation_nightly.csv'), index=False)

    def ba(x):
        x = x.dropna(); return pd.Series(dict(n=len(x), bias=x.mean(), median=x.median(), sd=x.std(), loa_low=x.mean() - 1.96 * x.std(),
                                              loa_high=x.mean() + 1.96 * x.std(), within30=(x.abs() <= 30).mean()))
    rows = []
    for meth, d in V[V.found == True].groupby('method'):
        rows.append(ba(d.onset_err_min).rename(lambda k: k).to_frame().T.assign(method=meth, measure='sleep onset'))
        rows.append(ba(d.wake_err_min).to_frame().T.assign(method=meth, measure='final awakening (all nights)'))
        rows.append(ba(d[d.ends_in_wake].wake_err_min).to_frame().T.assign(method=meth, measure='final awakening (Dreem ends in wake)'))
        if meth == 'model':
            for est, refc, lab_ in (('est_tst', 'ref_tst', 'TST'), ('est_deep', 'ref_n3', 'deep/N3 min'), ('est_rem', 'ref_rem', 'REM min')):
                rows.append(ba(d[est] - d[refc]).to_frame().T.assign(method=meth, measure=lab_ + ' (est - Dreem)'))
                rows.append(pd.DataFrame([dict(method=meth, measure=lab_ + ' r(est, Dreem)', n=int(d[[est, refc]].dropna().shape[0]),
                                               bias=float(d[est].corr(d[refc])))]))
    S = pd.concat(rows, ignore_index=True)[['method', 'measure', 'n', 'bias', 'median', 'sd', 'loa_low', 'loa_high', 'within30']]
    S.to_csv(os.path.join(R, 'spt_validation_summary.csv'), index=False)
    found = V.groupby('method').found.mean()
    G = V[(V.method == 'model') & V.found].assign(tst_err=lambda d: d.est_tst - d.ref_tst).groupby('genotype').agg(
        nights=('tst_err', 'size'), tst_bias=('tst_err', 'mean'), tst_sd=('tst_err', 'std'), onset_bias=('onset_err_min', 'mean'))
    G.to_csv(os.path.join(R, 'spt_validation_by_genotype.csv'))
    json.dump(dict(epoch_metrics=EP, spt_found_fraction=found.round(3).to_dict(), nights_scored=int(len(SP))), open(os.path.join(R, 'run_info.json'), 'w'), indent=1)
    pd.set_option('display.width', 220)
    print(json.dumps(EP)); print('SPT found on Dreem nights:', found.round(3).to_dict(), '| nights with model SPT:', len(SP))
    print(S.round(2).to_string(index=False)); print(G.round(1).to_string())
