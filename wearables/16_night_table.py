"""Night-level outcome table for the copy-number models, plus a large-error audit.

Rows = Dreem nights (sleep-architecture core QC pass, >= 80 % wrist coverage) - the 'paired' set -
and, separately, all wrist nights (retrained model 14) with >= 80 % of the sleep window worn.
Outcomes
  dreem_tst                Dreem TST (min)
  wrist_tst                wrist TST inside the wrist sleep window (retrained model)
  sol_headband             Dreem sleep onset - headband start (sleep-paper definition)
  sol_lightsout            Dreem sleep onset - lights-out estimate, where
                           lights-out = later of (mat in-bed, headband start)  [judgment call; see report]
  sol_matinbed             Dreem sleep onset - mat in-bed
  wrist_sol_lightsout      wrist sleep-window onset - mat in-bed (all mat nights, no headband needed)
  log versions: log(1 + SOL), as in the sleep paper.
Covariates: age_final, sex_final, site (redcap_data_access_group) from the 2026-10-12 manifest via linkage.
Usage: python -I 16_night_table.py <phase0_dir> <aim1_dir>
"""
import os, sys
import numpy as np, pandas as pd

P0, A1 = sys.argv[1:3]
R = os.path.join(A1, 'cn_models'); os.makedirs(R, exist_ok=True)
L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
cov = L[['empatica_id', 'base_id', 'guid', 'genotype', 'locus', 'copy_number', 'age_final', 'sex_final', 'redcap_data_access_group',
         'rarecnv_family_id', 'canonical_bp', 'link_status', 'visit']]
TR = pd.read_csv(os.path.join(A1, 'triangulation', 'triangulation_nightly.csv'))
for c in ['rec_start', 'rec_end', 'dreem_onset', 'dreem_final_wake', 'mat_in_bed', 'mat_out_bed', 'mat_final_sleep']:
    TR[c] = pd.to_datetime(TR[c], utc=True)
VV = pd.read_csv(os.path.join(A1, 'retrain_mat', 'spt_validation_nightly.csv'))
SP = pd.read_csv(os.path.join(A1, 'retrain_mat', 'spt_nightly.csv'))
for c in ['spt_on', 'spt_off']:
    SP[c] = pd.to_datetime(SP[c], utc=True)

# ---- paired (Dreem) nights ------------------------------------------------------------------------
N = TR[['recording_id', 'empatica_id', 'rec_start', 'rec_end', 'dreem_onset', 'dreem_tst', 'mat_in_bed', 'has_mat', 'post_class']].merge(
    VV[['recording_id', 'est_tst', 'end_err_triang', 'onset_err', 'ref_extended']], on='recording_id', how='left')
N = N.rename(columns={'est_tst': 'wrist_tst'})
# persistent sleep onset (sleep-architecture definition): start of the first 10-min window (20 epochs)
# that is >= 90 % headband-scored sleep (N1/N2/N3/REM; MT and wake count as not sleep)
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'), columns=['recording_id', 'epoch', 'ts_utc', 'stage'])
D = D[D.recording_id.isin(N.recording_id)].sort_values(['recording_id', 'epoch'])
def persist_onset(d, k=20):
    sl = d.stage.isin(['N1', 'N2', 'N3', 'R']).values.astype(float)
    if len(sl) < k:
        return pd.NaT
    run = np.convolve(sl, np.ones(k), 'valid')            # run[i] = sleep epochs in window starting at i
    ok = np.flatnonzero(run >= 0.9 * k)
    return d.ts_utc.iloc[ok[0]] if len(ok) else pd.NaT
po = D.groupby('recording_id').apply(persist_onset, include_groups=False).rename('dreem_persist_onset')
N = N.merge(po, left_on='recording_id', right_index=True, how='left')
N['dreem_persist_onset'] = pd.to_datetime(N.dreem_persist_onset, utc=True)
mins = lambda a, b: (a - b).dt.total_seconds() / 60
N['lightsout'] = N[['mat_in_bed', 'rec_start']].max(axis=1).where(N.mat_in_bed.notna())
# mat in-bed after Dreem sleep onset is impossible as lights-out (wrong session / mat detection delay):
# fall back to headband start for those nights, and record it
N['lightsout_fallback'] = N.lightsout.notna() & (N.lightsout > N.dreem_onset)
N.loc[N.lightsout_fallback, 'lightsout'] = N.loc[N.lightsout_fallback, 'rec_start']
N['sol_headband'] = mins(N.dreem_onset, N.rec_start)
N['solp_headband'] = mins(N.dreem_persist_onset, N.rec_start)
N['solp_lightsout'] = mins(N.dreem_persist_onset, N.lightsout)
N['solp_matinbed'] = mins(N.dreem_persist_onset, N.mat_in_bed)
for c in ['solp_headband', 'solp_lightsout', 'solp_matinbed']:
    N.loc[N[c] < 0, c] = np.nan
    N[c + '_log'] = np.log1p(N[c])
N['sol_lightsout'] = mins(N.dreem_onset, N.lightsout)
N['sol_matinbed'] = mins(N.dreem_onset, N.mat_in_bed)
for c in ['sol_headband', 'sol_lightsout', 'sol_matinbed']:
    N.loc[N[c] < 0, c] = np.nan if c == 'sol_matinbed' else N.loc[N[c] < 0, c]    # mat in-bed after onset: mis-matched session
    N[c + '_log'] = np.log1p(N[c].clip(lower=0))
N = N.merge(cov, on='empatica_id', how='left')
N.to_csv(os.path.join(R, 'nights_paired.csv'), index=False)

# ---- all wrist nights -------------------------------------------------------------------------------
rows = []
for pid, g in SP.groupby('empatica_id'):
    f = pd.read_parquet(os.path.join(A1, 'retrain_mat', f'scored_{pid}.parquet'))
    t = f.ts_utc.values
    for r in g.itertuples():
        m = (t >= np.datetime64(r.spt_on.tz_convert(None))) & (t < np.datetime64(r.spt_off.tz_convert(None)))
        rows.append(dict(empatica_id=pid, night=r.night, spt_on=r.spt_on, spt_off=r.spt_off, wrist_tst=r.tst_min,
                         spt_min=(r.spt_off - r.spt_on).total_seconds() / 60, worn_frac=float(f.worn.values[m].mean()) if m.any() else 0.0))
A = pd.DataFrame(rows)
WS = pd.read_parquet(os.path.join(P0, 'withings_sessions.parquet'))
A['study_id'] = A.empatica_id.str.replace('R2', '', regex=False)
mib = []
for r in A.itertuples():
    s = WS[(WS.study_id == r.study_id) & (WS.start_utc < r.spt_on + pd.Timedelta('1h')) & (WS.end_utc > r.spt_on)]
    mib.append(s.start_utc.max() if len(s) else pd.NaT)
A['mat_in_bed'] = pd.to_datetime(pd.Series(mib, index=A.index), utc=True)
A['wrist_sol_lightsout'] = mins(A.spt_on, A.mat_in_bed)
A.loc[A.wrist_sol_lightsout < 0, 'wrist_sol_lightsout'] = np.nan
A['wrist_sol_lightsout_log'] = np.log1p(A.wrist_sol_lightsout)
A = A.merge(cov, on='empatica_id', how='left')
A['valid'] = (A.worn_frac >= 0.8) & (A.spt_min >= 180)
A.to_csv(os.path.join(R, 'nights_wrist_all.csv'), index=False)

# ---- large-error audit (paired nights) ---------------------------------------------------------------
E = N.dropna(subset=['wrist_tst']).copy()
E['tst_err'] = E.wrist_tst - E.dreem_tst
E['big'] = E.tst_err.abs() > 120
def cause(r):
    if r.ref_extended is True or r.post_class == 'headband stopped, still asleep in bed':
        return 'headband truncated (mat: still asleep) - Dreem low, wrist right'
    if pd.notna(r.end_err_triang) and r.end_err_triang < -120:
        return 'wrist window ended early'
    if pd.notna(r.end_err_triang) and r.end_err_triang > 120:
        return 'wrist window ended late'
    if pd.notna(r.onset_err) and abs(r.onset_err) > 120:
        return 'wrist onset far from Dreem onset'
    return 'within-window misclassification (wake scored sleep or vice versa)'
E.loc[E.big, 'cause'] = E[E.big].apply(cause, axis=1)
au = E[E.big].groupby('cause').agg(nights=('tst_err', 'size'), median_err=('tst_err', 'median')).reset_index()
au.to_csv(os.path.join(R, 'large_error_audit.csv'), index=False)
pd.set_option('display.width', 200)
print('paired nights', len(N), '| with lights-out estimate', int(N.lightsout.notna().sum()), '| fallback to headband start', int(N.lightsout_fallback.sum()),
      '| negative SOL remaining:', {c: int((N[c] < 0).sum()) for c in ['sol_headband', 'sol_lightsout', 'sol_matinbed']})
print('persistent-sleep latency (min) quantiles: headband', N.solp_headband.quantile([.25, .5, .75]).round(1).tolist(),
      '| lights-out', N.solp_lightsout.quantile([.25, .5, .75]).round(1).tolist(), '| mat in-bed', N.solp_matinbed.quantile([.25, .5, .75]).round(1).tolist(),
      '| n', {c: int(N[c].notna().sum()) for c in ['solp_headband', 'solp_lightsout', 'solp_matinbed']},
      '| no persistent sleep:', int(N.dreem_persist_onset.isna().sum()))
print('first-epoch vs persistent onset gap (min) quantiles', mins(N.dreem_persist_onset, N.dreem_onset).quantile([.5, .75, .9, .99]).round(1).tolist())
print('SOL quantiles (min): headband', N.sol_headband.quantile([.25, .5, .75]).round(1).tolist(),
      '| lights-out', N.sol_lightsout.quantile([.25, .5, .75]).round(1).tolist(), '| mat in-bed', N.sol_matinbed.quantile([.25, .5, .75]).round(1).tolist())
print('all wrist nights', len(A), '| valid', int(A.valid.sum()), '| participants', A[A.valid].empatica_id.nunique(),
      '| valid with mat lights-out', int((A.valid & A.wrist_sol_lightsout.notna()).sum()))
print('large errors (|wrist - Dreem TST| > 120 min):', int(E.big.sum()), 'of', len(E)); print(au.round(1).to_string(index=False))
