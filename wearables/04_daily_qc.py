"""Phase 0.4: wear detection and valid-day QC.

Wear is assessed per 5-s epoch from two independent sources:
  * accelerometer non-wear (GGIR-style, van Hees 2011/2013): 60-min windows stepped every
    15 min; a window is non-wear when the pooled SD is < 13 mg on at least 2 of 3 axes.
    Pooled window SD is reconstructed from the 5-s epoch means and SDs
    (var = mean(within-epoch var) + var(epoch means)). GGIR also uses a value-range
    criterion, which needs the raw signal; it is not applied here.
    An epoch is non-wear if ANY window covering it is non-wear.
  * vendor wear flag: Empatica wearing_detection_percentage < 50 or a 'device_not_worn'
    missing reason in that minute.
  worn = data present AND acc-worn AND vendor-worn. Both components are kept so the
  combination rule can be changed in sensitivity analyses.
Days are local calendar days (midnight-midnight, participant time zone, DST-aware).
A valid day has >= H worn hours, H in {10, 16, 20}; primary H = 16.
A participant passes with >= 3 valid days (primary); the >= 1 weekend-day requirement
is reported as a sensitivity flag.

Usage: python -I 04_daily_qc.py <phase0_dir> [n_workers]
Writes <phase0>/qc/daily_qc.parquet, participant_qc.csv, qc_summary.csv, qc_venn_counts.csv
"""
import os, sys, glob, json
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

P0 = sys.argv[1]
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
QC = os.path.join(P0, 'qc'); os.makedirs(QC, exist_ok=True)
SD_TH, WIN, STEP = 0.013, 720, 180          # 13 mg; 60 min and 15 min in 5-s epochs
HOURS = (10, 16, 20)


def one(pid):
    f5 = os.path.join(P0, 'epoch5s', f'{pid}.parquet')
    fm = os.path.join(P0, 'minute', f'{pid}.parquet')
    if not os.path.exists(f5):
        return None
    e = pd.read_parquet(f5)
    m = pd.read_parquet(fm, columns=['ts_utc', 'wearing_detection_percentage', 'wear_missing_reason', 'tz']) if os.path.exists(fm) else None
    tz = m.tz.iloc[0] if m is not None and len(m) else 'UTC'
    # regular 5-s grid over the recording span; missing epochs = no data
    grid = pd.date_range(e.ts_utc.min().floor('D'), e.ts_utc.max().ceil('D'), freq='5s', tz='UTC', inclusive='left')
    e = e.set_index('ts_utc').reindex(grid)
    has = e.enmo.notna().values
    # pooled-window SD per axis via cumulative sums (NaN epochs contribute nothing)
    nonwear_win = np.zeros(len(e), bool)
    k = np.where(has, 1.0, 0.0)
    cs_n = np.r_[0, np.cumsum(k)]
    low_axes = np.zeros(len(e), int)
    starts = np.arange(0, max(len(e) - WIN + 1, 1), STEP)
    win_low = np.zeros(len(starts), int)
    for ax in 'xyz':
        mu = np.nan_to_num(e[ax].values); sd = np.nan_to_num(e['sd' + ax].values)
        c1 = np.r_[0, np.cumsum(mu * k)]; c2 = np.r_[0, np.cumsum(mu ** 2 * k)]; cv = np.r_[0, np.cumsum(sd ** 2 * k)]
        n = cs_n[starts + WIN] - cs_n[starts]
        with np.errstate(invalid='ignore', divide='ignore'):
            mean = (c1[starts + WIN] - c1[starts]) / n
            var = (cv[starts + WIN] - cv[starts]) / n + (c2[starts + WIN] - c2[starts]) / n - mean ** 2
        win_low += (np.sqrt(np.maximum(var, 0)) < SD_TH) & (n >= 0.5 * WIN)
    bad = starts[win_low >= 2]
    if len(bad):
        d = np.zeros(len(e) + 1, int)
        np.add.at(d, bad, 1); np.add.at(d, np.minimum(bad + WIN, len(e)), -1)
        nonwear_win = np.cumsum(d)[:-1] > 0
    acc_worn = has & ~nonwear_win
    # vendor wear flag, minute -> 5-s
    if m is not None and len(m):
        v = m.set_index('ts_utc')
        vw = (v.wearing_detection_percentage.fillna(0) >= 50) & ~v.wear_missing_reason.fillna('').str.contains('not_worn')
        vend = vw.reindex(grid.floor('min')).fillna(False).values.astype(bool)
        vend_known = v.wearing_detection_percentage.reindex(grid.floor('min')).notna().values | \
                     v.wear_missing_reason.reindex(grid.floor('min')).notna().values
    else:
        vend = np.ones(len(e), bool); vend_known = np.zeros(len(e), bool)
    worn = acc_worn & vend
    loc = grid.tz_convert(tz)
    df = pd.DataFrame({'date': loc.date, 'has': has, 'acc_worn': acc_worn, 'vend_worn': vend & vend_known, 'worn': worn,
                       'enmo': np.where(worn, e.enmo.values, np.nan)})
    g = df.groupby('date')
    day = pd.DataFrame({'data_h': g.has.sum() * 5 / 3600, 'acc_worn_h': g.acc_worn.sum() * 5 / 3600,
                        'vendor_worn_h': g.vend_worn.sum() * 5 / 3600, 'worn_h': g.worn.sum() * 5 / 3600,
                        'enmo_mg_worn': g.enmo.mean() * 1000}).reset_index()
    day['empatica_id'] = pid
    day['weekday'] = pd.to_datetime(day.date).dt.dayofweek
    day['weekend'] = day.weekday >= 5
    for h in HOURS:
        day[f'valid{h}'] = day.worn_h >= h
    return day


if __name__ == '__main__':
    L = pd.read_csv(os.path.join(P0, 'linkage.csv'))
    ids = sorted(os.path.basename(f)[:-8] for f in glob.glob(os.path.join(P0, 'epoch5s', '*.parquet')))
    with ProcessPoolExecutor(NW) as ex:
        days = pd.concat([d for d in ex.map(one, ids) if d is not None], ignore_index=True)
    days = days[days.data_h > 0].reset_index(drop=True)          # drop empty calendar days between recordings
    days.to_parquet(os.path.join(QC, 'daily_qc.parquet'), index=False)

    calib = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob(os.path.join(P0, 'calib', '*.json'))}
    pq = days.groupby('empatica_id').agg(qc_days=('date', 'size'), data_h=('data_h', 'sum'), worn_h=('worn_h', 'sum'),
                                         **{f'valid_days_{h}': (f'valid{h}', 'sum') for h in HOURS},
                                         valid_weekend_16=('valid16', lambda s: int((s & days.loc[s.index, 'weekend']).sum()))).reset_index()
    pq['calibrated'] = pq.empatica_id.map(lambda i: bool(calib.get(i, {}).get('converged', False))).astype(bool)
    pq = L.merge(pq, on='empatica_id', how='left')
    for h in HOURS:
        pq[f'pass_{h}h'] = pq[f'valid_days_{h}'].fillna(0) >= 3
    pq['pass_16h_weekend'] = pq.pass_16h & (pq.valid_weekend_16.fillna(0) >= 1)
    pq['in_analysis_groups'] = (pq.link_status == 'linked') & (pq.visit == 1)
    pq.to_csv(os.path.join(QC, 'participant_qc.csv'), index=False)

    real = pq[~pq.test_account]
    steps = [('Empatica IDs (excl. test accounts)', real),
             ('visit 1 (R2 repeats set aside)', real[real.visit == 1]),
             ('linked to a five-group genotype', real[real.in_analysis_groups]),
             ('raw accelerometer processed', real[real.in_analysis_groups & real.qc_days.notna()]),
             ('>= 3 valid days at 16 h (PRIMARY)', real[real.in_analysis_groups & real.pass_16h])]
    s = pd.DataFrame([(k, len(v), *[int((v.genotype == g).sum()) for g in ['16pDel', '16pDup', '22qDel', '22qDup', 'Control']])
                      for k, v in steps], columns=['step', 'n', '16pDel', '16pDup', '22qDel', '22qDup', 'Control'])
    a = real[real.in_analysis_groups]
    sens = pd.DataFrame([(f'>= 3 valid days at {h} h', int(a[f'pass_{h}h'].sum())) for h in HOURS] +
                        [('>= 3 valid days at 16 h incl. >= 1 weekend day', int(a.pass_16h_weekend.sum())),
                         ('calibration did not converge', int((a.calibrated == False).sum()))], columns=['criterion', 'n_pass'])
    s.to_csv(os.path.join(QC, 'qc_summary.csv'), index=False); sens.to_csv(os.path.join(QC, 'qc_sensitivity.csv'), index=False)
    # Venn of the three participant-level exclusion reasons among linked participants
    r = pd.DataFrame({'few_valid_days': ~a.pass_16h, 'no_weekend_day': a.pass_16h & ~a.pass_16h_weekend, 'uncalibrated': a.calibrated == False})
    r.groupby(list(r.columns)).size().reset_index(name='n').to_csv(os.path.join(QC, 'qc_venn_counts.csv'), index=False)
    print(s.to_string(index=False)); print(); print(sens.to_string(index=False))
    dd = days.merge(L[['empatica_id', 'genotype']], on='empatica_id')
    print('\nworn hours per day (all days with any data): median', round(days.worn_h.median(), 1),
          '| acc-worn but vendor-not-worn hours/day median', round((days.acc_worn_h - days.worn_h).median(), 2))
    print('valid16 days per participant by genotype (median):')
    print(dd.groupby(['genotype', 'empatica_id']).valid16.sum().groupby('genotype').median().to_string())
