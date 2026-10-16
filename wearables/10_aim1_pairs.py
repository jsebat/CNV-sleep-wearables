"""Aim 1, step 1: one row per Dreem 30-s epoch, with time-matched wrist features.

For every Dreem epoch [t, t+30 s) of every recording whose participant has Empatica data:
  accelerometer (5-s epochs inside the window): mean/max ENMO, mean anglez, max |d anglez|
  between consecutive 5-s epochs, n 5-s epochs present
  pulse/thermal (30-s epoch whose start is nearest t, within 15 s): HR, RMSSD, n beats, temp, EDA
  vendor (minute containing t): Empatica sleep_detection_stage, wearing %, activity counts VM
  wear (5-s epochs): fraction worn (accelerometer non-wear rule as in 04_daily_qc is NOT
  recomputed here; vendor wear % is carried and the Dreem window itself certifies the
  headband was on).
Context features (rolling, centred, over the same participant-night) are added in the model
step, not here, so this table stays a pure join.

Linking Dreem -> Empatica: Dreem file ID == Empatica study ID (R2 visits are separate
Empatica IDs; a Dreem recording is matched to whichever visit's data cover it).

Usage: python -I 10_aim1_pairs.py <phase0_dir> <out_dir>
"""
import os, sys
import numpy as np, pandas as pd

P0, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
ids = {r.empatica_id for r in L.itertuples() if not r.test_account}
rows = []
for rid, d in D.groupby('recording_id', sort=False):
    cand = [i for i in ids if i == d.dreem_id.iloc[0] or i == d.dreem_id.iloc[0] + 'R2']
    t0 = d.ts_utc.values.astype('datetime64[ns]')
    best = None
    for pid in cand:
        f5 = os.path.join(P0, 'epoch5s', f'{pid}.parquet')
        if not os.path.exists(f5):
            continue
        e = pd.read_parquet(f5, columns=['ts_utc', 'enmo', 'anglez'])
        te = e.ts_utc.values.astype('datetime64[ns]')
        i0 = np.searchsorted(te, t0); i1 = np.searchsorted(te, t0 + np.timedelta64(30, 's'))
        cover = (i1 > i0).mean()
        if best is None or cover > best[0]:
            best = (cover, pid, e, i0, i1)
    if best is None or best[0] == 0:
        continue
    cover, pid, e, i0, i1 = best
    en, an = e.enmo.values, e.anglez.values
    f = pd.DataFrame({'n5': i1 - i0})
    f['enmo_mean'] = [en[a:b].mean() if b > a else np.nan for a, b in zip(i0, i1)]
    f['enmo_max'] = [en[a:b].max() if b > a else np.nan for a, b in zip(i0, i1)]
    f['anglez_mean'] = [an[a:b].mean() if b > a else np.nan for a, b in zip(i0, i1)]
    f['dangle_max'] = [np.abs(np.diff(an[max(a - 1, 0):b])).max() if b - a > 0 and b > 1 else np.nan for a, b in zip(i0, i1)]
    x = d.reset_index(drop=True).assign(empatica_id=pid, **{c: f[c].values for c in f.columns})
    # 30-s physiology: nearest epoch start within 15 s
    f30 = os.path.join(P0, 'epoch30s', f'{pid}.parquet')
    if os.path.exists(f30):
        p = pd.read_parquet(f30).sort_values('ts_utc')
        p['ts_utc'] = p.ts_utc.astype('datetime64[ns, UTC]')
        x = pd.merge_asof(x.sort_values('ts_utc').assign(ts_utc=lambda z: z.ts_utc.astype('datetime64[ns, UTC]')), p,
                          on='ts_utc', direction='nearest', tolerance=pd.Timedelta('15s'))
    fm = os.path.join(P0, 'minute', f'{pid}.parquet')
    if os.path.exists(fm):
        m = pd.read_parquet(fm, columns=['ts_utc', 'sleep_detection_stage', 'wearing_detection_percentage', 'vector_magnitude'])
        m['ts_utc'] = m.ts_utc.astype('datetime64[ns, UTC]')
        x = pd.merge_asof(x.sort_values('ts_utc'), m.sort_values('ts_utc'), on='ts_utc', direction='backward',
                          tolerance=pd.Timedelta('60s'))
    x['coverage'] = cover
    rows.append(x)
P = pd.concat(rows, ignore_index=True)
P.to_parquet(os.path.join(OUT, 'aim1_epoch_pairs.parquet'), index=False)
r = P.groupby('recording_id').agg(id=('empatica_id', 'first'), cover=('coverage', 'first'), qc=('qc_pass_core', 'first'))
print('recordings with any wrist overlap', len(r), '| participants', r.id.nunique())
print('coverage quantiles', r.cover.quantile([.1, .25, .5]).round(3).tolist())
ok = r[(r.cover >= 0.8) & (r.qc == True)]
print('QC-pass recordings with >= 80 % wrist coverage:', len(ok), 'participants', ok.id.nunique())
print('epoch-level missingness:', P[['enmo_mean', 'hr_bpm', 'rmssd_ms', 'temp_c', 'sleep_detection_stage']].isna().mean().round(3).to_dict())
print(pd.crosstab(P.stage, P.sleep_detection_stage, dropna=False).to_string())
