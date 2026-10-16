"""Aim 4 second panel - EDA (skin conductance) from the raw 4 Hz WRIST BAND signal, one participant per call.

Sano, Picard & Stickgold 2014 (Int J Psychophysiol 94:382; doi 10.1016/j.ijpsycho.2014.09.011):
  low-pass 0.4 Hz (FIR, 32nd order); a PEAK is a rise whose first derivative exceeds 0.01 uS/s, peaks >= 1 s apart;
  a PEAK EPOCH is a 30-s epoch with >= 1 peak; a STORM is >= 2 consecutive peak epochs.
Per valid wrist night (sleep period from the retrained model, 20_phenotypes / nights_wrist_all), over sleep epochs:
  eda_peaks_per_h_sleep, eda_peak_epoch_frac (share of sleep epochs that are peak epochs), eda_storms_per_night,
  eda_storm_min (minutes in storms), eda_scl_sleep (median skin conductance level in sleep, uS)
Day (10:00-20:00, worn, awake): eda_scl_day (median), eda_peaks_per_h_day.
Usage: python -I 28_eda_raw.py <analysis_root> <empatica_root> <empatica_id>   -> aim4/phase2_eda/<ID>.csv
"""
import os, sys, glob
import numpy as np, pandas as pd, fastavro
from scipy.signal import firwin, filtfilt

ROOT, ER, PID = sys.argv[1:4]
P0, A1, A4 = (os.path.join(ROOT, d) for d in ('phase0', 'aim1', 'aim4'))
OUT = os.path.join(A4, 'phase2_eda'); os.makedirs(OUT, exist_ok=True)
fs = 4.0
T, V = [], []
for f in sorted(glob.glob(os.path.join(ER, '1', '1', 'participant_data', '*', f'{PID}-*', 'raw_data', 'v6', '*.avro'))):
    with open(f, 'rb') as fh:
        for rec in fastavro.reader(fh):
            e = rec['rawData']['eda']
            if len(e['values']) and e['samplingFrequency'] > 0:
                T.append(e['timestampStart'] / 1e6 + np.arange(len(e['values'])) / e['samplingFrequency']); V.append(np.asarray(e['values'], float))
if not T:
    sys.exit(0)
t = np.concatenate(T); v = np.concatenate(V); o = np.argsort(t); t, v = t[o], v[o]
keep = np.r_[True, np.diff(t) > 0]; t, v = t[keep], v[keep]
# filter within contiguous segments (gaps > 2 s split the signal)
seg = np.r_[0, np.flatnonzero(np.diff(t) > 2.0) + 1, len(t)]
b = firwin(33, 0.4, fs=fs)
vf = np.full_like(v, np.nan)
for a, z in zip(seg[:-1], seg[1:]):
    if z - a > 3 * len(b):
        vf[a:z] = filtfilt(b, [1.0], v[a:z])
d = np.r_[np.nan, np.diff(vf) * fs]                       # uS / s
cand = np.flatnonzero((d > 0.01) & ~np.isnan(d))
# a peak = onset of a rise above threshold; successive candidates within 1 s belong to the same peak
peaks = cand[np.r_[True, np.diff(t[cand]) >= 1.0]] if len(cand) else cand
pt = pd.to_datetime(t[peaks], unit='s', utc=True)
s30 = pd.Series(1, index=pt).resample('30s').sum()
scl = pd.Series(vf, index=pd.to_datetime(t, unit='s', utc=True)).resample('30s').median()
f30 = pd.read_parquet(os.path.join(P0, 'features30s', f'{PID}.parquet'), columns=['ts_utc', 'ts_local', 'worn'])
sc = pd.read_parquet(os.path.join(A1, 'retrain_mat', f'scored_{PID}.parquet'), columns=['ts_utc', 'p_sleep'])
f30 = f30.merge(sc, on='ts_utc', how='left').set_index('ts_utc')
f30['peaks'] = s30.reindex(f30.index).fillna(0); f30['scl'] = scl.reindex(f30.index)
f30['has_eda'] = scl.reindex(f30.index).notna()
N = pd.read_csv(os.path.join(A1, 'cn_models', 'nights_wrist_all.csv'))
N = N[(N.empatica_id == PID) & N.valid.astype(str).isin(['True', '1'])]
rows = []
for r in N.itertuples():
    w = f30.loc[pd.Timestamp(r.spt_on):pd.Timestamp(r.spt_off)]
    w = w[w.worn & (w.p_sleep >= 0.5) & w.has_eda]
    if len(w) < 120:
        continue
    pe = (w.peaks >= 1).values
    dd = np.diff(np.r_[0, pe.astype(int), 0]); st, en = np.flatnonzero(dd == 1), np.flatnonzero(dd == -1)
    storms = (en - st) >= 2
    rows.append(dict(empatica_id=PID, night=r.night, eda_peaks_per_h_sleep=w.peaks.sum() / (len(w) / 120),
                     eda_peak_epoch_frac=pe.mean(), eda_storms_per_night=int(storms.sum()),
                     eda_storm_min=float(((en - st)[storms]).sum() / 2), eda_scl_sleep=w.scl.median()))
day = f30[f30.worn & (f30.p_sleep < 0.5) & f30.has_eda & (f30.ts_local.dt.hour >= 10) & (f30.ts_local.dt.hour < 20)]
out = pd.DataFrame(rows)
out['eda_scl_day'] = day.scl.median() if len(day) else np.nan
out['eda_peaks_per_h_day'] = day.peaks.sum() / (len(day) / 120) if len(day) else np.nan
out.to_csv(os.path.join(OUT, f'{PID}.csv'), index=False)
print(PID, 'nights', len(out), 'peaks', len(peaks))
