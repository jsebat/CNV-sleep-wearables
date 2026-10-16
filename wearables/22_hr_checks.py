"""Aim 4, heart-rate confounder checks - part 1: per-night signal-quality and still-epoch measures (WRIST BAND),
plus heart-related history from REDCap.

Per valid wrist night (same sleep periods as 20_phenotypes), over epochs the retrained model scores as sleep:
  hr_coverage      share of sleep epochs with a heart-rate value (pulse beats detected and clean)
  beats_per_epoch  mean clean beats per 30-s sleep epoch
  enmo_sleep_mg    mean wrist movement during sleep (mg)
  still_frac       share of sleep epochs with max 5-s arm-angle change <= 5 deg and ENMO < 10 mg
  hr_still_mean    mean HR in still sleep epochs only;  rmssd_still  median RMSSD in still sleep epochs
REDCap (Aug 2026): med_birthdefects_3 (congenital malformation of the heart), med_birthdefects_4 (circulatory),
  med_cardiovascular (any cardiovascular problem) - 1 = yes.
Usage: python -I 22_hr_checks.py <analysis_root> <redcap_csv> [n_workers]
"""
import os, sys, csv
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

ROOT, RC = sys.argv[1:3]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
P0, A1, A4 = (os.path.join(ROOT, d) for d in ('phase0', 'aim1', 'aim4'))
NIGHTS = pd.read_csv(os.path.join(A1, 'cn_models', 'nights_wrist_all.csv'))
NIGHTS = NIGHTS[NIGHTS.valid.astype(str).isin(['True', '1'])]
for c in ['spt_on', 'spt_off']:
    NIGHTS[c] = pd.to_datetime(NIGHTS[c], utc=True)


def one(pid):
    fp = os.path.join(P0, 'features30s', f'{pid}.parquet'); fs = os.path.join(A1, 'retrain_mat', f'scored_{pid}.parquet')
    if not (os.path.exists(fp) and os.path.exists(fs)):
        return []
    f = pd.read_parquet(fp, columns=['ts_utc', 'worn', 'enmo_mean', 'dangle_max', 'hr_bpm', 'rmssd_ms', 'n_beats'])
    f = f.merge(pd.read_parquet(fs, columns=['ts_utc', 'p_sleep']), on='ts_utc', how='left')
    out = []
    for r in NIGHTS[NIGHTS.empatica_id == pid].itertuples():
        seg = f[(f.ts_utc >= r.spt_on) & (f.ts_utc < r.spt_off) & f.worn & (f.p_sleep >= 0.5)]
        if len(seg) < 120:
            continue
        still = (seg.dangle_max <= 5) & (seg.enmo_mean < 0.010)
        st = seg[still]
        out.append(dict(empatica_id=pid, night=r.night, hr_coverage=seg.hr_bpm.notna().mean(), beats_per_epoch=seg.n_beats.mean(),
                        enmo_sleep_mg=1000 * seg.enmo_mean.mean(), still_frac=still.mean(),
                        hr_still_mean=st.hr_bpm.mean() if len(st) >= 60 else np.nan,
                        rmssd_still=st.rmssd_ms.median() if len(st) >= 60 else np.nan))
    return out


if __name__ == '__main__':
    ids = sorted(NIGHTS.empatica_id.unique())
    with ProcessPoolExecutor(NW) as ex:
        N = pd.DataFrame([x for xs in ex.map(one, ids) for x in xs])
    N.to_csv(os.path.join(A4, 'hr_checks_night.csv'), index=False)
    L = pd.read_csv(os.path.join(P0, 'linkage.csv')); ids_r = set(L.rarecnv_id.dropna())
    want = ['rarecnv_id', 'redcap_event_name', 'med_birthdefects_3', 'med_birthdefects_4', 'med_cardiovascular']
    rows = []
    with open(RC, newline='', encoding='utf-8', errors='replace') as fh:
        rd = csv.reader(fh); h = next(rd); ix = {c: h.index(c) for c in want if c in h}
        print('REDCap heart fields found:', [c for c in want if c in ix])
        for x in rd:
            if x[ix['rarecnv_id']] in ids_r:
                rows.append({c: x[i] for c, i in ix.items()})
    R = pd.DataFrame(rows).replace({'': np.nan})
    H = R.groupby('rarecnv_id').agg({c: (lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan) for c in ix if c not in ('rarecnv_id', 'redcap_event_name')})
    H = L[['empatica_id', 'rarecnv_id', 'genotype']].merge(H, left_on='rarecnv_id', right_index=True, how='left')
    H.to_csv(os.path.join(A4, 'heart_history.csv'), index=False)
    pd.set_option('display.width', 200)
    print('nights', len(N), '| participants', N.empatica_id.nunique())
    M = N.merge(L[['empatica_id', 'genotype']], on='empatica_id')
    print(M.groupby('genotype')[['hr_coverage', 'beats_per_epoch', 'enmo_sleep_mg', 'still_frac']].median().round(3).to_string())
    for c in [c for c in ix if c.startswith('med_')]:
        print(c, H.groupby('genotype')[c].apply(lambda s: f"{(s == '1').sum()}/{s.notna().sum()}").to_dict())
