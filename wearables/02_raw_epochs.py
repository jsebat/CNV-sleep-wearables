"""Phase 0c: raw EmbracePlus Avro -> epoch-level signals, one participant per call.

Accelerometer (64 Hz, ADC counts -> g via the file's imuParams):
  * auto-calibration (van Hees et al. 2014 sphere fit on still 10-s windows; offset + scale per axis)
  * 5-s epochs: ENMO = mean(max(|a| - 1, 0)) [van Hees 2013], anglez (GGIR definition: angle of the
    z axis from 5-s rolling median of x,y,z), per-axis SD (for non-wear), mean x/y/z
Pulse (systolic peak times -> inter-beat intervals; kept if 300-2000 ms and within 20 % of the
  previous IBI), 30-s epochs: mean HR, RMSSD, n clean beats.
Temperature (1 Hz) and EDA (4 Hz): 30-s means.
All timestamps are UTC; local time is joined later from the minute tables (same tz).

Usage: python -I 02_raw_epochs.py <empatica_root> <out_dir> <empatica_id>
Writes <out>/epoch5s/<ID>.parquet, <out>/epoch30s/<ID>.parquet, <out>/calib/<ID>.json
"""
import os, sys, glob, json
import numpy as np, pandas as pd
import fastavro

ROOT, OUT, PID = sys.argv[1:4]
for d in ('epoch5s', 'epoch30s', 'calib'):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)
files = sorted(glob.glob(os.path.join(ROOT, '1', '1', 'participant_data', '*', f'{PID}-*', 'raw_data', 'v6', '*.avro')))


def adc_to_g(v, p):
    return (np.asarray(v, float) - p['digitalMin']) * (p['physicalMax'] - p['physicalMin']) / \
           (p['digitalMax'] - p['digitalMin']) + p['physicalMin']


acc, ibi, tmp, eda = [], [], [], []
meta = dict(n_files=len(files), fw=set(), algo=set(), skipped_fs0=dict(acc=0, temperature=0, eda=0), skipped_fs0_samples=dict(acc=0, temperature=0, eda=0))
for f in files:
    with open(f, 'rb') as fh:
        for rec in fastavro.reader(fh):
            meta['fw'].add('%d.%d.%d' % tuple(rec['fwVersion'].values()))
            meta['algo'].add('%d.%d.%d' % tuple(rec['algoVersion'].values()))
            r = rec['rawData']
            a = r['accelerometer']
            if len(a['x']) and not (a['samplingFrequency'] > 0):
                meta['skipped_fs0']['acc'] += 1; meta['skipped_fs0_samples']['acc'] += len(a['x'])
            elif len(a['x']):
                t = a['timestampStart'] / 1e6 + np.arange(len(a['x'])) / a['samplingFrequency']
                p = a['imuParams']
                acc.append(np.column_stack([t, adc_to_g(a['x'], p), adc_to_g(a['y'], p), adc_to_g(a['z'], p)]))
            pk = np.asarray(r['systolicPeaks']['peaksTimeNanos'], float) / 1e9
            if len(pk) > 1:
                ibi.append(np.column_stack([pk[1:], np.diff(pk) * 1000]))
            for sig, store in (('temperature', tmp), ('eda', eda)):
                s = r[sig]
                if len(s['values']) and not (s['samplingFrequency'] > 0):
                    meta['skipped_fs0'][sig] += 1; meta['skipped_fs0_samples'][sig] += len(s['values'])
                elif len(s['values']):
                    t = s['timestampStart'] / 1e6 + np.arange(len(s['values'])) / s['samplingFrequency']
                    store.append(np.column_stack([t, s['values']]))
if not acc:
    print(PID, 'no accelerometer data'); sys.exit(0)
A = np.concatenate(acc); A = A[np.isfinite(A).all(axis=1)]; A = A[np.argsort(A[:, 0], kind='stable')]
A = A[np.concatenate([[True], np.diff(A[:, 0]) > 0])]           # drop overlapping samples between files
fs = 64.0

# ---- auto-calibration (sphere fit on still windows) ----------------------------------
w = int(10 * fs)
n = len(A) // w
blk = A[:n * w, 1:].reshape(n, w, 3)
still = (blk.std(axis=1) < 0.013).all(axis=1)
M = blk.mean(axis=1)[still]
off, scl = np.zeros(3), np.ones(3)
calib = dict(n_still_windows=int(still.sum()), converged=False)
if len(M) >= 50 and (M.max(0) - M.min(0)).min() > 0.3:      # need spread on every axis to identify scale
    for _ in range(100):
        C = (M + off) * scl
        tgt = C / np.linalg.norm(C, axis=1, keepdims=True)
        for k in range(3):                                  # per-axis linear regression tgt ~ a + b*raw
            X = np.column_stack([np.ones(len(M)), M[:, k]])
            b = np.linalg.lstsq(X, tgt[:, k], rcond=None)[0]
            scl[k], off[k] = b[1], b[0] / b[1]
        err = np.abs(np.linalg.norm((M + off) * scl, axis=1) - 1).mean()
    calib.update(converged=bool(err < 0.01), error_after=float(err),
                 error_before=float(np.abs(np.linalg.norm(M, axis=1) - 1).mean()))
    if err < 0.01:
        A[:, 1:] = (A[:, 1:] + off) * scl
calib.update(offset=off.tolist(), scale=scl.tolist(), fw=sorted(meta['fw']), algo=sorted(meta['algo']),
             n_files=meta['n_files'], n_samples=int(len(A)), skipped_fs0_segments=meta['skipped_fs0'],
             skipped_fs0_samples=meta['skipped_fs0_samples'])
json.dump(calib, open(os.path.join(OUT, 'calib', f'{PID}.json'), 'w'), indent=1)

# ---- 5-s epochs -----------------------------------------------------------------------
df = pd.DataFrame(A, columns=['t', 'x', 'y', 'z'])
roll = df[['x', 'y', 'z']].rolling(int(5 * fs), center=True, min_periods=1).median()
df['anglez'] = np.degrees(np.arctan(roll.z / np.sqrt(roll.x ** 2 + roll.y ** 2)))
df['enmo'] = np.maximum(np.sqrt(df.x ** 2 + df.y ** 2 + df.z ** 2) - 1, 0)
df['ep'] = (df.t // 5).astype('int64')
g = df.groupby('ep')
e5 = pd.DataFrame({'enmo': g.enmo.mean(), 'anglez': g.anglez.mean(), 'x': g.x.mean(), 'y': g.y.mean(), 'z': g.z.mean(),
                   'sdx': g.x.std(), 'sdy': g.y.std(), 'sdz': g.z.std(), 'n': g.t.size()})
e5 = e5[e5.n >= 0.8 * 5 * fs]
e5.index = pd.to_datetime(e5.index * 5, unit='s', utc=True); e5.index.name = 'ts_utc'
e5.astype({c: 'float32' for c in e5.columns if c != 'n'}).reset_index().to_parquet(os.path.join(OUT, 'epoch5s', f'{PID}.parquet'), index=False)

# ---- 30-s epochs: pulse, temperature, EDA ---------------------------------------------
out30 = []
if ibi:
    I = np.concatenate(ibi); I = I[np.argsort(I[:, 0])]
    ok = (I[:, 1] > 300) & (I[:, 1] < 2000)
    ok[1:] &= np.abs(np.diff(I[:, 1])) / I[:-1, 1] < 0.2
    b = pd.DataFrame({'t': I[:, 0], 'ibi': I[:, 1], 'ok': ok})
    b['sd'] = np.r_[np.nan, np.diff(b.ibi)] ** 2
    b.loc[~(b.ok & b.ok.shift(1, fill_value=False)), 'sd'] = np.nan
    b = b[b.ok]; b['ep'] = (b.t // 30).astype('int64')
    gb = b.groupby('ep')
    out30.append(pd.DataFrame({'hr_bpm': 60000 / gb.ibi.mean(), 'rmssd_ms': np.sqrt(gb.sd.mean()), 'n_beats': gb.ibi.size()}))
for name, store in (('temp_c', tmp), ('eda_us', eda)):
    if store:
        S = np.concatenate(store)
        out30.append(pd.DataFrame({name: S[:, 1]}, index=(S[:, 0] // 30).astype('int64')).groupby(level=0).mean())
if out30:
    e30 = pd.concat(out30, axis=1)
    e30.index = pd.to_datetime(e30.index * 30, unit='s', utc=True); e30.index.name = 'ts_utc'
    e30.astype('float32').reset_index().to_parquet(os.path.join(OUT, 'epoch30s', f'{PID}.parquet'), index=False)
print(PID, 'files', len(files), 'epochs5s', len(e5), 'calib', calib.get('converged'), round(calib.get('error_before', np.nan), 4),
      '->', round(calib.get('error_after', np.nan), 4))
