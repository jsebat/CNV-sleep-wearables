"""Core HEADBAND (Dreem EEG) sleep phenotypes per participant, for the master phenotype table.

Nights: recordings passing the sleep-architecture core QC (782 recordings). Per-recording values are taken from the
sleep-architecture table (SleepArchitectureResults per_recording_metrics_with_qc.csv) so they match that analysis
exactly: TST, TRT, SE, WASO, SOL (from headband start), SOL_persist, min/pct N1 N2 N3 REM, DI = N3/(N1+N2+REM),
DI_over_NREM = N3/(N1+N2+N3), DI_over_TST = N3/TST, MT fraction.
Added here (from the 30-s hypnograms):
  rem_latency_min      sleep onset (first N1/N2/N3/REM epoch) to the first REM epoch (AASM convention)
  n_awakenings         wake runs >= 2 epochs (1 min) between sleep onset and the last sleep epoch
  sol_lightsout_min    sleep onset - lights-out (later of mat in-bed and headband start; 16_night_table), where available
Per participant: mean over QC-pass nights (n_hb_nights recorded). Participants are keyed by the Dreem/Empatica study ID.
Usage: python -I 29_headband_core.py <analysis_root> -> aim4/headband_core_participant.csv, headband_core_night.csv
"""
import os, sys
import numpy as np, pandas as pd

ROOT = sys.argv[1]
P0, A1, A4 = (os.path.join(ROOT, d) for d in ('phase0', 'aim1', 'aim4'))
Q = pd.read_csv(os.path.join(ROOT, 'inputs', 'per_recording_metrics_with_qc.csv'))
Q = Q[Q.qc_pass_core == True]
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'), columns=['recording_id', 'dreem_id', 'epoch', 'stage'])
D = D[D.recording_id.isin(Q.recording_id)].sort_values(['recording_id', 'epoch'])
SL = {'N1', 'N2', 'N3', 'R'}
rows = []
for rid, d in D.groupby('recording_id', sort=False):
    s = d.stage.values; sl = np.isin(s, list(SL))
    if not sl.any():
        continue
    on = np.argmax(sl); last = len(sl) - 1 - np.argmax(sl[::-1])
    rem = np.flatnonzero(s == 'R')
    w = (s[on:last + 1] == 'W').astype(int); dd = np.diff(np.r_[0, w, 0]); st, en = np.flatnonzero(dd == 1), np.flatnonzero(dd == -1)
    rows.append(dict(recording_id=rid, dreem_id=d.dreem_id.iloc[0], rem_latency_min=(rem[rem >= on][0] - on) / 2 if (rem >= on).any() else np.nan,
                     n_awakenings=int(((en - st) >= 2).sum())))
X = pd.DataFrame(rows)
N = Q.merge(X, on='recording_id')
NP = pd.read_csv(os.path.join(A1, 'cn_models', 'nights_paired.csv'), usecols=['recording_id', 'sol_lightsout'])
N = N.merge(NP, on='recording_id', how='left')
keep = ['TST', 'TRT', 'SE', 'WASO', 'SOL', 'SOL_persist', 'min_N1', 'min_N2', 'min_N3', 'min_R', 'pct_N1', 'pct_N2', 'pct_N3', 'pct_R',
        'DI', 'DI_over_NREM', 'DI_over_TST', 'frac_MT', 'rem_latency_min', 'n_awakenings', 'sol_lightsout']
N[['recording_id', 'dreem_id'] + keep].to_csv(os.path.join(A4, 'headband_core_night.csv'), index=False)
P = N.groupby('dreem_id')[keep].mean()
P['n_hb_nights'] = N.groupby('dreem_id').size()
P['n_hb_nights_lightsout'] = N.groupby('dreem_id').sol_lightsout.count()
P = P.rename(columns={c: 'hb_' + c for c in keep}).reset_index().rename(columns={'dreem_id': 'study_id'})
P.to_csv(os.path.join(A4, 'headband_core_participant.csv'), index=False)
pd.set_option('display.width', 200)
print('headband participants', len(P), '| nights', len(N))
print(P.drop(columns='study_id').describe(percentiles=[.5]).T[['count', 'mean', '50%', 'min', 'max']].round(2).to_string())
