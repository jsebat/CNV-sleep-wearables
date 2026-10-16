"""Aim 1b: triangulate Dreem (EEG headband), Withings (under-mattress) and the wrist model.

Per Dreem night (QC-pass, >= 80 % wrist coverage) with an overlapping Withings bed session:
 1. clock alignment: minute-level sleep/wake agreement Dreem vs Withings at lags -60..+60 min;
    the best lag should be ~0 if both clocks are right.
 2. arbitration of the post-headband interval (Dreem recording end .. wrist SPT end):
      in_bed_frac  = share of minutes the mat reports a state (person on the mat)
      asleep_frac  = share of minutes with mat state > 0
    classes: 'headband stopped, still asleep in bed' (in bed >= 50 % and asleep >= 50 % of in-bed),
             'in bed awake', 'out of bed', 'no mat data'.
 3. headband truncation: mat out-of-bed time minus Dreem recording end, mat final sleep minus Dreem
    final sleep; by genotype.
 4. three-way night agreement: TST, time in bed, sleep onset latency from lights-out (mat in-bed time).
Usage: python -I 13_triangulate.py <phase0_dir> <aim1_dir>
"""
import os, sys
import numpy as np, pandas as pd

P0, A1 = sys.argv[1:3]
R = os.path.join(A1, 'triangulation'); os.makedirs(R, exist_ok=True)
L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
WS = pd.read_parquet(os.path.join(P0, 'withings_sessions.parquet'))
WT = pd.read_parquet(os.path.join(P0, 'withings_states.parquet'))
V = pd.read_csv(os.path.join(A1, 'deployable', 'spt_validation_nightly.csv'))
V = V[(V.method == 'model') & (V.found == True)]
SP = pd.read_csv(os.path.join(A1, 'deployable', 'spt_model_nightly.csv'), parse_dates=['spt_on', 'spt_off'])
SLEEP = ['N1', 'N2', 'N3', 'R']
out, lagrows = [], []
for r in V.itertuples():
    d = D[D.recording_id == r.recording_id]
    sid = d.dreem_id.iloc[0]
    rec_start, rec_end = d.ts_utc.min(), d.ts_utc.max() + pd.Timedelta('30s')
    sl = d.stage.isin(SLEEP)
    on, fw = d.ts_utc[sl].min(), d.ts_utc[sl].max() + pd.Timedelta('30s')
    s = WS[(WS.study_id == sid) & (WS.start_utc < rec_end) & (WS.end_utc > rec_start)]
    rec = dict(recording_id=r.recording_id, empatica_id=r.empatica_id, study_id=sid, rec_start=rec_start, rec_end=rec_end,
               dreem_onset=on, dreem_final_wake=fw, dreem_tst=sl.sum() / 2, wrist_tst=r.est_tst, has_mat=len(s) > 0)
    c = SP[(SP.empatica_id == r.empatica_id) & (SP.spt_on < fw) & (SP.spt_off > on)]
    if len(c):
        c = c.iloc[(np.minimum(c.spt_off, fw) - np.maximum(c.spt_on, on)).values.argmax()]
        rec.update(wrist_on=c.spt_on, wrist_off=c.spt_off)
    if len(s):
        s = s.iloc[(np.minimum(s.end_utc, rec_end) - np.maximum(s.start_utc, rec_start)).values.argmax()]
        st = WT[(WT.study_id == sid) & (WT.ts_utc >= s.start_utc - pd.Timedelta('2h')) & (WT.ts_utc <= s.end_utc + pd.Timedelta('2h'))].set_index('ts_utc').state
        ms = st[st > 0]
        rec.update(mat_in_bed=s.start_utc, mat_out_bed=s.end_utc, mat_tst=(s.light_s + s.deep_s + s.rem_s) / 60,
                   mat_to_sleep=s.to_sleep_s / 60, mat_final_sleep=ms.index.max() + pd.Timedelta('1min') if len(ms) else pd.NaT,
                   mat_onset=ms.index.min() if len(ms) else pd.NaT)
        # 1. alignment: Dreem minute sleep vs mat minute sleep, lags
        dm = d.assign(m=d.ts_utc.dt.floor('min')).groupby('m').stage.agg(lambda x: x.isin(SLEEP).mean() >= 0.5)
        best = (np.nan, -1)
        for lag in range(-60, 61, 2):
            mt = st.copy(); mt.index = mt.index + pd.Timedelta(minutes=lag)
            j = dm.to_frame('dreem').join((mt > 0).rename('mat'), how='inner')
            if len(j) >= 120:
                a = (j.dreem == j.mat).mean()
                if a > best[1]:
                    best = (lag, a)
                if lag == 0:
                    rec['agree_lag0'] = a
        rec.update(best_lag=best[0], agree_best=best[1])
        # 2. arbitration of the post-headband interval
        if 'wrist_off' in rec and rec['wrist_off'] > rec_end + pd.Timedelta('10min'):
            mins = pd.date_range(rec_end.ceil('min'), rec['wrist_off'].floor('min'), freq='min')
            present = st.reindex(mins)
            inbed = present.notna().mean(); asleep = (present > 0).sum() / max(present.notna().sum(), 1)
            cls = ('headband stopped, still asleep in bed' if inbed >= 0.5 and asleep >= 0.5 else
                   'in bed awake' if inbed >= 0.5 else 'out of bed')
            rec.update(post_min=len(mins), post_in_bed_frac=inbed, post_asleep_frac=asleep, post_class=cls)
        elif 'wrist_off' in rec:
            rec.update(post_min=0, post_class='no post-headband interval')
    out.append(rec)
T = pd.DataFrame(out).merge(L[['empatica_id', 'genotype', 'age_final']], on='empatica_id', how='left')
for a, b, n in (('mat_out_bed', 'rec_end', 'mat_out_minus_headband_end'), ('mat_final_sleep', 'dreem_final_wake', 'mat_final_sleep_minus_dreem'),
                ('mat_in_bed', 'rec_start', 'mat_in_minus_headband_start'), ('wrist_off', 'mat_final_sleep', 'wrist_off_minus_mat_final_sleep'),
                ('dreem_onset', 'mat_in_bed', 'dreem_sol_from_mat_lights_out')):
    T[n] = (T[a] - T[b]).dt.total_seconds() / 60
T.to_csv(os.path.join(R, 'triangulation_nightly.csv'), index=False)
M = T[T.has_mat]
pd.set_option('display.width', 220)
print('Dreem nights', len(T), '| with overlapping mat session', len(M), '| participants', M.empatica_id.nunique())
print('\n1. clock alignment: best lag (min) quantiles', M.best_lag.quantile([.1, .25, .5, .75, .9]).tolist(),
      '| agreement at lag 0 median %.3f, at best lag %.3f' % (M.agree_lag0.median(), M.agree_best.median()))
print('\n2. post-headband interval (nights where wrist SPT ends > 10 min after headband):')
print(M.post_class.value_counts().to_string())
pm = M[M.post_min > 0]
print('minutes in each class:', pm.groupby('post_class').post_min.sum().round(0).to_dict())
q = lambda s: s.quantile([.1, .25, .5, .75, .9]).round(1).tolist()
print('\n3. truncation: mat out-of-bed minus headband end (min)', q(M.mat_out_minus_headband_end))
print('   mat final sleep minus Dreem final sleep (min)', q(M.mat_final_sleep_minus_dreem))
print('   wrist SPT end minus mat final sleep (min)', q(M.wrist_off_minus_mat_final_sleep))
print('   mat in-bed minus headband start (min)', q(M.mat_in_minus_headband_start))
g = M.groupby('genotype').agg(nights=('recording_id', 'size'), trunc_median=('mat_final_sleep_minus_dreem', 'median'),
                              trunc_gt60=('mat_final_sleep_minus_dreem', lambda x: (x > 60).mean()),
                              still_asleep_class=('post_class', lambda x: (x == 'headband stopped, still asleep in bed').mean()))
print(g.round(2).to_string())
print('\n4. TST (min) medians: Dreem %.0f, mat %.0f, wrist %.0f' % (M.dreem_tst.median(), M.mat_tst.median(), M.wrist_tst.median()))
for a, b in (('mat_tst', 'dreem_tst'), ('wrist_tst', 'dreem_tst'), ('wrist_tst', 'mat_tst')):
    x = (M[a] - M[b]).dropna(); print(f'   {a} - {b}: bias {x.mean():.1f}, sd {x.std():.1f}, r {M[a].corr(M[b]):.2f}, n {len(x)}')
print('   Dreem SOL from mat lights-out (min) quantiles', q(M.dreem_sol_from_mat_lights_out), '| mat own latency', q(M.mat_to_sleep))
