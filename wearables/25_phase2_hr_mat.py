"""Aim 4, second panel - block 1: WRIST BAND heart-rate rhythm and MAT (Withings) breathing / heart rate.

WRIST BAND (per person, valid days/nights from 20_phenotypes / 04_daily_qc):
  hr_day_still      mean HR in awake, still, worn daytime epochs (10:00-20:00 local; P(sleep) < 0.5;
                    max 5-s arm-angle change <= 5 deg; ENMO < 10 mg) - a resting-HR proxy [definition provisional
                    until the citation check returns]
  hr_dip_pct        100 * (hr_day_still - mean sleep HR) / hr_day_still   (night-time dip)
  hr_cos_mesor / hr_cos_amp / hr_cos_acrophase_h   single 24-h cosinor on hourly mean worn HR over valid days
MAT (per mat sleep session matched to the participant; session = Withings bed session; minutes with a sleep state):
  mat_hr_sleep      mean mat HR over minutes the mat scores as asleep (README: raw_bed_hr, bpm)
  mat_rr_sleep      mean respiratory rate over asleep minutes (README: raw_bed_respiratory-rate, breaths/min)
  mat_snore_frac    share of in-bed minutes with snoring = 1 (README: raw_bed_snoring, 0/1)
  mat_breath_score_pos_frac   EXPLORATORY - share of in-bed minutes with 'Breathing event probability' > 0;
                    undocumented vendor score, not a validated apnea index
Usage: python -I 25_phase2_hr_mat.py <analysis_root> <withings_root> [n_workers]
Writes aim4/phase2_hr_mat_participant.csv, aim4/phase2_mat_session.csv
"""
import os, sys, glob, json, re
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

ROOT, WR = sys.argv[1:3]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
P0, A1, A4 = (os.path.join(ROOT, d) for d in ('phase0', 'aim1', 'aim4'))
DAYS = pd.read_parquet(os.path.join(P0, 'qc', 'daily_qc.parquet')); DAYS = DAYS[DAYS.valid16]
PN = pd.read_csv(os.path.join(A4, 'phenotypes_night.csv'))


def wrist(pid):
    fp = os.path.join(P0, 'features30s', f'{pid}.parquet'); fs = os.path.join(A1, 'retrain_mat', f'scored_{pid}.parquet')
    if not (os.path.exists(fp) and os.path.exists(fs)):
        return None
    f = pd.read_parquet(fp, columns=['ts_utc', 'ts_local', 'worn', 'enmo_mean', 'dangle_max', 'hr_bpm'])
    f = f.merge(pd.read_parquet(fs, columns=['ts_utc', 'p_sleep']), on='ts_utc', how='left')
    vd = set(DAYS[DAYS.empatica_id == pid].date.astype(str))
    f = f[f.worn & f.ts_local.dt.date.astype(str).isin(vd)]
    if len(vd) < 3 or f.hr_bpm.notna().sum() < 1000:
        return None
    h = f.ts_local.dt.hour
    still_awake_day = (h >= 10) & (h < 20) & (f.p_sleep < 0.5) & (f.dangle_max <= 5) & (f.enmo_mean < 0.010)
    day = f.loc[still_awake_day, 'hr_bpm'].mean()
    sl = PN[PN.empatica_id == pid].hr_sleep_mean.mean()
    hh = f.set_index('ts_local').hr_bpm.resample('1h').mean().dropna()
    t = (hh.index.hour + 0.5).values
    X = np.column_stack([np.ones(len(t)), np.cos(2 * np.pi * t / 24), np.sin(2 * np.pi * t / 24)])
    b = np.linalg.lstsq(X, hh.values, rcond=None)[0]
    return dict(empatica_id=pid, hr_day_still=day, n_day_still_epochs=int(still_awake_day.sum()),
                hr_dip_pct=100 * (day - sl) / day if np.isfinite(day) and np.isfinite(sl) else np.nan,
                hr_cos_mesor=b[0], hr_cos_amp=float(np.hypot(b[1], b[2])), hr_cos_acrophase_h=float((np.degrees(np.arctan2(b[2], b[1])) / 15) % 24))


IDRE = re.compile(r'(MONT|TORO|PENN|UCLA)(\d{4,5})'); MONTRE = re.compile(r'211-\d+[_-](3\d{3})')
def sid_of(path):
    m = IDRE.search(path)
    if m:
        return m.group(1) + m.group(2)
    m = MONTRE.search(path); return 'MONT' + m.group(1) if m else None


def expand(fn, sid):
    r = pd.read_csv(fn); rows = []
    for t0, du, va in zip(r.start, r.duration, r.value):
        t = pd.Timestamp(t0).tz_convert('UTC'); du = json.loads(du); va = json.loads(va)
        for a, v in zip(du, va):
            for k in range(max(int(a) // 60, 1)):
                rows.append((t + pd.Timedelta(minutes=k), v))
            t += pd.Timedelta(seconds=a)
    return pd.DataFrame(rows, columns=['ts_utc', 'value']).assign(study_id=sid)


def mat_signals():
    roots = [WR, os.path.join(P0, 'withings_extracted')]
    out = {}
    for name, key in (('raw_bed_hr.csv', 'hr'), ('raw_bed_respiratory-rate.csv', 'rr'), ('raw_bed_snoring.csv', 'snore'),
                      ('raw_bed_Breathing event probability.csv', 'breath')):
        parts = []
        for root in roots:
            for fn in glob.glob(os.path.join(root, '**', name), recursive=True):
                sid = sid_of(os.path.relpath(fn, root))
                if sid and os.path.getsize(fn) > 30:
                    parts.append(expand(fn, sid))
        out[key] = pd.concat(parts).drop_duplicates(['study_id', 'ts_utc']) if parts else pd.DataFrame(columns=['ts_utc', 'value', 'study_id'])
    return out


if __name__ == '__main__':
    L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
    ids = sorted(L[~L.test_account].empatica_id)
    with ProcessPoolExecutor(NW) as ex:
        Wr = pd.DataFrame([r for r in ex.map(wrist, ids) if r])
    S = mat_signals()
    WS = pd.read_parquet(os.path.join(P0, 'withings_sessions.parquet')); WT = pd.read_parquet(os.path.join(P0, 'withings_states.parquet'))
    rows = []
    for sid, ss in WS.groupby('study_id'):
        st = WT[WT.study_id == sid].set_index('ts_utc').state
        sig = {k: v[v.study_id == sid].set_index('ts_utc').value for k, v in S.items()}
        for s in ss.itertuples():
            inbed = st[(st.index >= s.start_utc) & (st.index < s.end_utc)]
            if len(inbed) < 180:
                continue
            asleep = inbed.index[inbed > 0]
            g = lambda k, idx: sig[k].reindex(idx)
            rows.append(dict(study_id=sid, start_utc=s.start_utc, minutes_in_bed=len(inbed), minutes_asleep=len(asleep),
                             mat_hr_sleep=g('hr', asleep).mean(), mat_rr_sleep=g('rr', asleep).mean(),
                             mat_snore_frac=(g('snore', inbed.index) == 1).mean() if g('snore', inbed.index).notna().any() else np.nan,
                             mat_breath_score_pos_frac=(g('breath', inbed.index) > 0).mean() if g('breath', inbed.index).notna().any() else np.nan))
    MS = pd.DataFrame(rows); MS.to_csv(os.path.join(A4, 'phase2_mat_session.csv'), index=False)
    MP = MS.groupby('study_id')[['mat_hr_sleep', 'mat_rr_sleep', 'mat_snore_frac', 'mat_breath_score_pos_frac']].mean().reset_index()
    MP['n_mat_sessions'] = MS.groupby('study_id').size().values
    link = L[['empatica_id', 'base_id']].rename(columns={'base_id': 'study_id'})
    P = Wr.merge(link, on='empatica_id', how='outer').merge(MP, on='study_id', how='left')
    P.to_csv(os.path.join(A4, 'phase2_hr_mat_participant.csv'), index=False)
    pd.set_option('display.width', 200)
    print('wrist participants', len(Wr), '| mat sessions', len(MS), '| mat participants', MS.study_id.nunique())
    print(P.drop(columns=['empatica_id', 'study_id']).describe(percentiles=[.5]).T[['count', 'mean', '50%', 'min', 'max']].round(2).to_string())
    # cross-device check: mat HR vs wrist HR in sleep, same people
    x = P.merge(PN.groupby('empatica_id').hr_sleep_mean.mean().reset_index(), on='empatica_id')
    print('r(mat sleep HR, wrist sleep HR) =', round(x.mat_hr_sleep.corr(x.hr_sleep_mean), 3), 'n =', int(x[['mat_hr_sleep', 'hr_sleep_mean']].dropna().shape[0]))
