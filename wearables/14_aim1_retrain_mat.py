"""Aim 1, step 4: retrain the wrist sleep model with mat-certified out-of-bed wake examples.

Why: the deployable model (12) was trained only on Dreem-recorded epochs, i.e. almost all in bed,
and its sleep window ran a median ~1 h past the mat's out-of-bed time (13_triangulate).
Training data here = (a) Dreem-labelled epochs as in 12, plus (b) wrist epochs the Withings mat
shows as OUT OF BED, labelled wake (y4 = W): the 3 h after each mat out-of-bed time and the 2 h
before each in-bed time, excluding any minute the mat records a state (someone on the mat) and
any non-worn epoch. Same features (features30s), same XGBoost settings, 5 folds grouped by family.
Evaluation
  * Dreem-labelled epochs only (comparable with 12): accuracy, sensitivity, specificity, kappa, 4-class.
  * held-out out-of-bed epochs: share scored wake.
  * sleep window (same rule as 12) vs references: onset vs Dreem onset; end vs the triangulated
    final awakening = Dreem final sleep, replaced by the mat's final sleep when the mat shows the
    person still asleep in bed after the headband stopped (13 class 'headband stopped, still asleep').
Usage: python -I 14_aim1_retrain_mat.py <phase0_dir> <aim1_dir> [n_threads]
"""
import os, sys, json
import numpy as np, pandas as pd
import xgboost as xgb
from sklearn.model_selection import GroupKFold
from sklearn.metrics import cohen_kappa_score, f1_score

P0, A1 = sys.argv[1:3]
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 12
FT = os.path.join(P0, 'features30s')
R = os.path.join(A1, 'retrain_mat'); os.makedirs(R, exist_ok=True)
L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
feats = json.load(open(os.path.join(A1, 'deployable', 'features.json')))


def runs(b):
    d = np.diff(np.r_[0, b.astype(int), 0]); return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def longest_block(b, min_len, max_gap):
    s, e = runs(b)
    if not len(s):
        return None
    bl = [[s[0], e[0]]]
    for a, z in zip(s[1:], e[1:]):
        if a - bl[-1][1] < max_gap:
            bl[-1][1] = z
        else:
            bl.append([a, z])
    bl = [x for x in bl if x[1] - x[0] >= min_len]
    return max(bl, key=lambda x: x[1] - x[0]) if bl else None


# ---- (a) Dreem-labelled epochs (as in 12) ---------------------------------------------------------
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
pairs = pd.read_parquet(os.path.join(A1, 'aim1_epoch_pairs.parquet'), columns=['recording_id', 'empatica_id', 'coverage']).drop_duplicates('recording_id')
keep = pairs[pairs.coverage >= 0.8]
DL = D.merge(keep[['recording_id', 'empatica_id']], on='recording_id')
DL = DL[(DL.qc_pass_core == True) & ~DL.stage.isin(['MT', 'U'])].copy()
DL['ts_grid'] = DL.ts_utc.dt.round('30s')
# ---- (b) mat out-of-bed windows ---------------------------------------------------------------------
WS = pd.read_parquet(os.path.join(P0, 'withings_sessions.parquet'))
WT = pd.read_parquet(os.path.join(P0, 'withings_states.parquet'))
emp_ids = set(L.empatica_id)
parts = []
for pid in sorted(emp_ids):
    fp = os.path.join(FT, f'{pid}.parquet')
    if not os.path.exists(fp):
        continue
    sid = pid.replace('R2', '')
    f = pd.read_parquet(fp)
    lab = pd.Series(np.nan, index=f.index); src = pd.Series('', index=f.index); rid = pd.Series('', index=f.index)
    d = DL[DL.empatica_id == pid].drop_duplicates('ts_grid')
    if len(d):
        ix = f.ts_utc.isin(d.ts_grid)
        m = d.set_index('ts_grid')
        lab[ix] = f.loc[ix, 'ts_utc'].map(m.stage.map({'W': 0, 'N1': 1, 'N2': 1, 'N3': 2, 'R': 3}))
        rid[ix] = f.loc[ix, 'ts_utc'].map(m.recording_id); src[ix] = 'dreem'
    s = WS[WS.study_id == sid]
    if len(s):
        on_mat = set(WT[WT.study_id == sid].ts_utc.dt.floor('min'))
        win = np.zeros(len(f), bool); t = f.ts_utc.values
        for a, b in zip(s.end_utc, s.start_utc):
            win |= (t >= np.datetime64(a.tz_convert(None))) & (t < np.datetime64((a + pd.Timedelta('3h')).tz_convert(None)))
            win |= (t >= np.datetime64((b - pd.Timedelta('2h')).tz_convert(None))) & (t < np.datetime64(b.tz_convert(None)))
        win &= ~f.ts_utc.dt.floor('min').isin(on_mat).values & f.worn.values & lab.isna().values
        # never inside any mat session
        for a, b in zip(s.start_utc, s.end_utc):
            win &= ~((t >= np.datetime64(a.tz_convert(None))) & (t < np.datetime64(b.tz_convert(None))))
        lab[win] = 0; src[win] = 'mat_out_of_bed'
    sel = lab.notna() & f.worn
    if sel.any():
        parts.append(f.loc[sel].assign(y4=lab[sel].astype(int).values, src=src[sel].values, recording_id=rid[sel].values, empatica_id=pid))
T = pd.concat(parts, ignore_index=True)
T = T.merge(L[['empatica_id', 'base_id', 'rarecnv_family_id', 'genotype']], on='empatica_id', how='left')
T['family'] = T.rarecnv_family_id.fillna(T.base_id)
T['sleep'] = (T.y4 > 0).astype(int)
print('training epochs by source:', T.src.value_counts().to_dict(), '| participants', T.empatica_id.nunique(), '| families', T.family.nunique())

X = T[feats].astype('float32').values
prm = dict(n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, tree_method='hist', n_jobs=NT)
T['p_sleep'] = np.nan; T['pred4'] = np.nan; T['fold'] = -1
for k, (tr, te) in enumerate(GroupKFold(5).split(X, T.sleep, groups=T.family)):
    T.loc[te, 'p_sleep'] = xgb.XGBClassifier(**prm, eval_metric='logloss').fit(X[tr], T.sleep.values[tr]).predict_proba(X[te])[:, 1]
    T.loc[te, 'pred4'] = xgb.XGBClassifier(**prm, objective='multi:softprob', num_class=4).fit(X[tr], T.y4.values[tr]).predict(X[te])
    T.loc[te, 'fold'] = k
Dt = T[T.src == 'dreem']; yb = (Dt.p_sleep >= 0.5).astype(int)
f1 = f1_score(Dt.y4, Dt.pred4.astype(int), average=None, labels=[0, 1, 2, 3])
EP = dict(dreem_epochs=len(Dt), accuracy=float((yb == Dt.sleep).mean()), sensitivity=float(yb[Dt.sleep == 1].mean()),
          specificity=float(1 - yb[Dt.sleep == 0].mean()), kappa=float(cohen_kappa_score(Dt.sleep, yb)),
          acc4=float((Dt.y4 == Dt.pred4).mean()), kappa4=float(cohen_kappa_score(Dt.y4, Dt.pred4.astype(int))),
          F1_W_light_deep_REM=[round(float(v), 3) for v in f1],
          out_of_bed_epochs=int((T.src == 'mat_out_of_bed').sum()),
          out_of_bed_scored_wake=float((T[T.src == 'mat_out_of_bed'].p_sleep < 0.5).mean()))
final_b = xgb.XGBClassifier(**prm, eval_metric='logloss').fit(X, T.sleep.values)
final_4 = xgb.XGBClassifier(**prm, objective='multi:softprob', num_class=4).fit(X, T.y4.values)
final_b.save_model(os.path.join(R, 'model_binary.json')); final_4.save_model(os.path.join(R, 'model_4class.json'))

# ---- score all epochs (out-of-fold where labelled), sleep window per night -------------------------
oof = T.set_index(['empatica_id', 'ts_utc'])[['p_sleep', 'pred4']]
SP = []
for pid in sorted(emp_ids):
    fp = os.path.join(FT, f'{pid}.parquet')
    if not os.path.exists(fp):
        continue
    f = pd.read_parquet(fp)
    Xf = f[feats].astype('float32').values
    f['p_sleep'] = final_b.predict_proba(Xf)[:, 1]; f['pred4'] = final_4.predict(Xf)
    if pid in oof.index.get_level_values(0):
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
            SP.append(dict(empatica_id=pid, night=n, spt_on=seg.ts_utc.iloc[0], spt_off=seg.ts_utc.iloc[-1] + pd.Timedelta('30s'),
                           tst_min=float((seg.p_sleep >= 0.5).sum() / 2),
                           **{f'min_{k}': float((seg.pred4 == i).sum() / 2) for i, k in enumerate(['W', 'light', 'deep', 'REM'])}))
SP = pd.DataFrame(SP); SP.to_csv(os.path.join(R, 'spt_nightly.csv'), index=False)

# ---- validation against Dreem and the triangulated final awakening ----------------------------------
TR = pd.read_csv(os.path.join(A1, 'triangulation', 'triangulation_nightly.csv'),
                 parse_dates=['dreem_onset', 'dreem_final_wake', 'mat_final_sleep', 'mat_out_bed', 'rec_end'])
for c in ['dreem_onset', 'dreem_final_wake', 'mat_final_sleep', 'mat_out_bed', 'rec_end']:
    TR[c] = pd.to_datetime(TR[c], utc=True)
TR['ref_final_wake'] = TR.dreem_final_wake
still = TR.post_class == 'headband stopped, still asleep in bed'
TR.loc[still, 'ref_final_wake'] = TR.loc[still, 'mat_final_sleep']
TR['ref_extended'] = still
rows = []
for r in TR.itertuples():
    c = SP[(SP.empatica_id == r.empatica_id) & (SP.spt_on < r.ref_final_wake) & (SP.spt_off > r.dreem_onset)]
    if not len(c):
        rows.append(dict(recording_id=r.recording_id, found=False)); continue
    c = c.iloc[(np.minimum(c.spt_off, r.ref_final_wake) - np.maximum(c.spt_on, r.dreem_onset)).values.argmax()]
    rows.append(dict(recording_id=r.recording_id, empatica_id=r.empatica_id, genotype=r.genotype, found=True, has_mat=r.has_mat,
                     ref_extended=r.ref_extended, onset_err=(c.spt_on - r.dreem_onset).total_seconds() / 60,
                     end_err_dreem=(c.spt_off - r.dreem_final_wake).total_seconds() / 60,
                     end_err_triang=(c.spt_off - r.ref_final_wake).total_seconds() / 60,
                     end_minus_mat_out=(c.spt_off - r.mat_out_bed).total_seconds() / 60 if pd.notna(r.mat_out_bed) else np.nan,
                     est_tst=c.tst_min, dreem_tst=r.dreem_tst))
VV = pd.DataFrame(rows); VV.to_csv(os.path.join(R, 'spt_validation_nightly.csv'), index=False)
F = VV[VV.found == True]


def summ(x, name):
    x = x.dropna(); return dict(measure=name, n=len(x), bias=x.mean(), median=x.median(), sd=x.std(),
                                loa_low=x.mean() - 1.96 * x.std(), loa_high=x.mean() + 1.96 * x.std(), within30=(x.abs() <= 30).mean())
S = pd.DataFrame([summ(F.onset_err, 'sleep onset vs Dreem'),
                  summ(F.end_err_dreem, 'window end vs Dreem final sleep'),
                  summ(F.end_err_triang, 'window end vs triangulated final sleep'),
                  summ(F[F.has_mat == True].end_err_triang, 'window end vs triangulated (mat nights)'),
                  summ(F.end_minus_mat_out, 'window end minus mat out-of-bed'),
                  summ(F.est_tst - F.dreem_tst, 'TST est - Dreem (all)'),
                  summ((F.est_tst - F.dreem_tst)[F.ref_extended == False], 'TST est - Dreem (headband not truncated)')])
S.to_csv(os.path.join(R, 'spt_validation_summary.csv'), index=False)
G = F.assign(d=F.est_tst - F.dreem_tst).groupby('genotype').agg(nights=('d', 'size'), tst_bias=('d', 'mean'), tst_sd=('d', 'std'),
                                                                end_err=('end_err_triang', 'median'))
G.to_csv(os.path.join(R, 'spt_validation_by_genotype.csv'))
json.dump(dict(epoch=EP, found=float(VV.found.mean()), nights_scored=len(SP)), open(os.path.join(R, 'run_info.json'), 'w'), indent=1)
pd.set_option('display.width', 220)
print(json.dumps(EP)); print('SPT found on Dreem nights %.3f | nights scored %d' % (VV.found.mean(), len(SP)))
print(S.round(2).to_string(index=False)); print(G.round(1).to_string())
