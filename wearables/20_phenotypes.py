"""Aim 4, step 1: standard actigraphy phenotypes - WRIST BAND (Empatica EmbracePlus) only.

Inputs (phase0 / aim1): features30s/<ID> (30-s ENMO, anglez, HR, RMSSD, temp, worn, local time, night),
  retrain_mat/scored_<ID> (P(sleep) from the retrained wrist model), cn_models/nights_wrist_all.csv
  (wrist sleep period per night, valid flag), qc/daily_qc.parquet (valid days, >= 16 h worn).
Levels
  NIGHT (valid wrist nights: >= 80 % of the sleep period worn, period >= 3 h)
    spt_h, tst_h, sleep_eff (tst/spt), waso_min, n_wake_bouts (>= 2 consecutive wake epochs inside period),
    onset_clock_h, offset_clock_h, midpoint_clock_h (local hours since previous noon, 12 = midnight... see below),
    hr_sleep_min (5th percentile of 30-s HR in sleep epochs), hr_sleep_mean, rmssd_sleep (median), temp_sleep_mean
    clock times are expressed as hours after 18:00 local (e.g. 23:30 -> 5.5, 07:00 -> 13.0) to avoid the midnight wrap.
  DAY (valid calendar days, >= 16 h worn)
    enmo_mg (mean worn ENMO), sed_h / light_h / mvpa_min (Hildebrand 2014 adult non-dominant wrist cut-points
    44.8 / 100.6 mg; per-minute ENMO; SECONDARY - not calibrated for this device), m10_mg, l5_mg,
    m10_onset_h, l5_onset_h (clock hours), temp_amp (90th - 10th percentile of worn 30-s wrist temperature)
  PARTICIPANT (needs >= 3 valid days)
    IS, IV (Van Someren 1999 definitions on hourly mean ENMO over valid days), RA = (M10 - L5) / (M10 + L5) from the
    average 24-h profile, cosinor mesor / amplitude / acrophase (single 24-h cosinor on hourly log10 ENMO),
    SRI (Phillips 2017: 200 * P(same sleep/wake state 24 h apart) - 100, over worn epoch pairs, 30-s resolution),
    tst_sd_h, midpoint_sd_h (night-to-night SD), social_jetlag_h (|mean midpoint Fri/Sat nights - other nights|),
    plus participant means of every night and day measure.
Usage: python -I 20_phenotypes.py <phase0_dir> <aim1_dir> <out_dir> [n_workers]
"""
import os, sys, glob
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

P0, A1, OUT = sys.argv[1:4]
NW = int(sys.argv[4]) if len(sys.argv) > 4 else 8
os.makedirs(OUT, exist_ok=True)
NIGHTS = pd.read_csv(os.path.join(A1, 'cn_models', 'nights_wrist_all.csv'))
NIGHTS = NIGHTS[NIGHTS.valid.astype(str).isin(['True', '1'])]
for c in ['spt_on', 'spt_off']:
    NIGHTS[c] = pd.to_datetime(NIGHTS[c], utc=True)
DAYS = pd.read_parquet(os.path.join(P0, 'qc', 'daily_qc.parquet'))
DAYS = DAYS[DAYS.valid16]
SED, MOD = 44.8, 100.6           # mg, Hildebrand 2014 adult non-dominant wrist (ActiGraph)


def clock18(ts):                  # local timestamp -> hours after 18:00 (18:00 = 0, 06:00 = 12)
    h = ts.dt.hour + ts.dt.minute / 60
    return (h - 18) % 24


def one(pid):
    fp = os.path.join(P0, 'features30s', f'{pid}.parquet')
    fs = os.path.join(A1, 'retrain_mat', f'scored_{pid}.parquet')
    if not (os.path.exists(fp) and os.path.exists(fs)):
        return None
    f = pd.read_parquet(fp, columns=['ts_utc', 'ts_local', 'worn', 'enmo_mean', 'hr_bpm', 'rmssd_ms', 'temp_c'])
    s = pd.read_parquet(fs, columns=['ts_utc', 'p_sleep'])
    f = f.merge(s, on='ts_utc', how='left')
    f['sleep'] = (f.p_sleep >= 0.5).where(f.worn & f.p_sleep.notna())
    f['date'] = f.ts_local.dt.date
    # ---------------- nights ----------------
    nrows = []
    for r in NIGHTS[NIGHTS.empatica_id == pid].itertuples():
        seg = f[(f.ts_utc >= r.spt_on) & (f.ts_utc < r.spt_off)]
        if not len(seg):
            continue
        sl = seg.sleep.fillna(False).astype(bool).values
        wake = ~sl
        d = np.diff(np.r_[0, wake.astype(int), 0]); starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        bouts = int(((ends - starts) >= 2).sum())
        ss = seg[sl]
        on_l, off_l = seg.ts_local.iloc[0], seg.ts_local.iloc[-1]
        on18, off18 = clock18(pd.Series([on_l])).iloc[0], clock18(pd.Series([off_l])).iloc[0]
        if off18 < on18:
            off18 += 24
        nrows.append(dict(empatica_id=pid, night=r.night, weekday=pd.Timestamp(on_l - pd.Timedelta(hours=12)).dayofweek,
                          spt_h=len(seg) / 120, tst_h=sl.sum() / 120, sleep_eff=sl.mean(), waso_min=wake.sum() / 2,
                          n_wake_bouts=bouts, onset_clock_h=on18, offset_clock_h=off18, midpoint_clock_h=(on18 + off18) / 2,
                          hr_sleep_min=ss.hr_bpm.quantile(0.05), hr_sleep_mean=ss.hr_bpm.mean(),
                          rmssd_sleep=ss.rmssd_ms.median(), temp_sleep_mean=ss.temp_c.mean()))
    N = pd.DataFrame(nrows)
    # ---------------- days ----------------
    vd = set(DAYS[DAYS.empatica_id == pid].date.astype(str))
    w = f[f.worn & f.date.astype(str).isin(vd)].copy()
    drows = []
    m = w.set_index('ts_local').enmo_mean.resample('1min').mean() * 1000      # per-minute ENMO (mg), worn only
    for dte, mm in m.groupby(m.index.date):
        if str(dte) not in vd or mm.notna().sum() < 16 * 60:
            continue
        roll10 = mm.rolling(600, min_periods=480).mean(); roll5 = mm.rolling(300, min_periods=240).mean()
        wd = w[w.date == dte]
        drows.append(dict(empatica_id=pid, date=str(dte), weekday=pd.Timestamp(dte).dayofweek, enmo_mg=mm.mean(),
                          sed_h=(mm < SED).sum() / 60, light_h=((mm >= SED) & (mm < MOD)).sum() / 60, mvpa_min=float((mm >= MOD).sum()),
                          m10_mg=roll10.max(), l5_mg=roll5.min(),
                          m10_onset_h=((roll10.idxmax() - pd.Timedelta(minutes=599)).hour + (roll10.idxmax() - pd.Timedelta(minutes=599)).minute / 60) if roll10.notna().any() else np.nan,
                          l5_onset_h=((roll5.idxmin() - pd.Timedelta(minutes=299)).hour + (roll5.idxmin() - pd.Timedelta(minutes=299)).minute / 60) if roll5.notna().any() else np.nan,
                          temp_amp=wd.temp_c.quantile(0.9) - wd.temp_c.quantile(0.1)))
    D = pd.DataFrame(drows)
    # ---------------- participant ----------------
    P = dict(empatica_id=pid, n_valid_days=len(D), n_valid_nights=len(N))
    if len(D) >= 3:
        h = m.resample('1h').mean().dropna()
        h = h[pd.Index(h.index.date).astype(str).isin(set(D.date))]
        x = h.values; n = len(x); xm = x.mean()
        prof = h.groupby(h.index.hour).mean()
        P['IS'] = (n * ((prof - xm) ** 2).sum()) / (24 * ((x - xm) ** 2).sum())
        P['IV'] = (n * (np.diff(x) ** 2).sum()) / ((n - 1) * ((x - xm) ** 2).sum())
        prof_m = m[pd.Index(m.index.date).astype(str).isin(set(D.date))]
        pm = prof_m.groupby(prof_m.index.hour * 60 + prof_m.index.minute).mean()
        circ = pd.concat([pm, pm]).reset_index(drop=True)
        m10 = circ.rolling(600).mean().max(); l5 = circ.rolling(300).mean().min()
        P['RA'] = (m10 - l5) / (m10 + l5)
        y = np.log10(h.clip(lower=0.1)).values; t = (h.index.hour + 0.5).values
        X = np.column_stack([np.ones(len(t)), np.cos(2 * np.pi * t / 24), np.sin(2 * np.pi * t / 24)])
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        P['cos_mesor'] = b[0]; P['cos_amp'] = float(np.hypot(b[1], b[2])); P['cos_acrophase_h'] = float((np.degrees(np.arctan2(b[2], b[1])) / 15) % 24)
        for c in ['enmo_mg', 'sed_h', 'light_h', 'mvpa_min', 'm10_mg', 'l5_mg', 'm10_onset_h', 'l5_onset_h', 'temp_amp']:
            P[c] = D[c].mean()
    st = f.set_index('ts_utc').sleep
    lag = st.shift(freq='24h').reindex(st.index)
    ok = st.notna() & lag.notna()
    if ok.sum() > 2 * 24 * 120:
        P['SRI'] = 200 * (st[ok].astype(bool) == lag[ok].astype(bool)).mean() - 100
    if len(N) >= 3:
        for c in ['spt_h', 'tst_h', 'sleep_eff', 'waso_min', 'n_wake_bouts', 'onset_clock_h', 'offset_clock_h', 'midpoint_clock_h',
                  'hr_sleep_min', 'hr_sleep_mean', 'rmssd_sleep', 'temp_sleep_mean']:
            P[c] = N[c].mean()
        P['tst_sd_h'] = N.tst_h.std(); P['midpoint_sd_h'] = N.midpoint_clock_h.std()
        we = N.weekday.isin([4, 5])                                         # Friday and Saturday nights
        if we.any() and (~we).any():
            P['social_jetlag_h'] = abs(N[we].midpoint_clock_h.mean() - N[~we].midpoint_clock_h.mean())
    return P, N, D


if __name__ == '__main__':
    ids = sorted(os.path.basename(p)[:-8] for p in glob.glob(os.path.join(P0, 'features30s', '*.parquet')) if not p.endswith('_build_info.json'))
    with ProcessPoolExecutor(NW) as ex:
        res = [r for r in ex.map(one, ids) if r is not None]
    P = pd.DataFrame([r[0] for r in res]); N = pd.concat([r[1] for r in res if len(r[1])]); D = pd.concat([r[2] for r in res if len(r[2])])
    P.to_csv(os.path.join(OUT, 'phenotypes_participant.csv'), index=False)
    N.to_csv(os.path.join(OUT, 'phenotypes_night.csv'), index=False)
    D.to_csv(os.path.join(OUT, 'phenotypes_day.csv'), index=False)
    pd.set_option('display.width', 220)
    print('participants', len(P), '| with >= 3 valid days', int((P.n_valid_days >= 3).sum()), '| nights', len(N), '| days', len(D))
    print(P.drop(columns=['empatica_id']).describe(percentiles=[.5]).T[['count', 'mean', '50%', 'min', 'max']].round(2).to_string())
