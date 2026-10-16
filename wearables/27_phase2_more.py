"""Aim 4 second panel - blocks 2-5, WRIST BAND only (kept separate from the core phenotypes).
Valid days = >= 16 h worn (04_daily_qc); per-person values = mean over valid days / nights; >= 3 valid days required.

ACTIVITY STRUCTURE (5-s ENMO, worn epochs of valid days)
  ig            intensity gradient (Rowlands 2018, MSSE 50:1323): 25-mg ENMO bins 0-4000 mg, minutes per bin,
                ln(minutes) ~ ln(bin mid-point) slope; empty bins dropped (GGIR g.intensitygradient).
  m60/m30/m15/m5_mg  MX: acceleration above which the most active X minutes of the day are accrued (Rowlands 2019).
  astp, satp    1 / mean active-bout and 1 / mean sedentary-bout duration (min^-1) (Schrack 2019, J Gerontol A 74:560),
                1-min ENMO; ACTIVE = >= 44.8 mg (Hildebrand 2014 light-activity cut-point) - our adaptation of the
                count threshold to ENMO.
NAPS
  daytime_sib_min    per day, outside the wrist sleep period: minutes in sustained-inactivity bouts (van Hees 2015:
                     >= 5 min with no 5-s arm-angle change > 5 deg); daytime_sib30_min = bouts >= 30 min only (the
                     30-min bout length of Dashti 2021; their inactivity threshold is not stated in the verified text).
CHRONOTYPE (Roenneberg 2019, Biology 8:54)
  msf_h         mid-sleep on free days (Fri/Sat nights) = SOf + SDf/2;  msfsc_h = MSF if SDf <= SDw, else SOf + SDweek/2,
                SDweek = (5 SDw + 2 SDf)/7 [weighting is common usage, not written out in the 2019 review].
                SO = wrist sleep-period onset, SD = sleep-period length; clock = hours after 18:00.
STEPS (Empatica per-minute step counts, worn minutes of valid days)
  steps_day; peak30_cadence = mean of the 30 highest step-count minutes of each day (Tudor-Locke 2012, JPAH 9:1125).
WRIST TEMPERATURE (Martinez-Nicolas 2019, Front Physiol 10:1396; 2013 PLoS One 8:e61142)
  10-min mean temperature profile averaged over valid days: temp_m5 (warmest 5 h), temp_l10 (coolest 10 h),
  temp_m5_onset_h, temp_l10_onset_h, temp_ra = (M5 - L10)/(M5 + L10); cosinor temp_cos_mesor, temp_cos_acrophase_h.
Usage: python -I 27_phase2_more.py <analysis_root> [n_workers] -> aim4/phase2_blocks2to5_participant.csv
"""
import os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

ROOT = sys.argv[1]; NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
P0, A1, A4 = (os.path.join(ROOT, d) for d in ('phase0', 'aim1', 'aim4'))
DAYS = pd.read_parquet(os.path.join(P0, 'qc', 'daily_qc.parquet')); DAYS = DAYS[DAYS.valid16]
NIGHTS = pd.read_csv(os.path.join(A1, 'cn_models', 'nights_wrist_all.csv'))
NIGHTS = NIGHTS[NIGHTS.valid.astype(str).isin(['True', '1'])]
for c in ['spt_on', 'spt_off']:
    NIGHTS[c] = pd.to_datetime(NIGHTS[c], utc=True)
PN = pd.read_csv(os.path.join(A4, 'phenotypes_night.csv'))


def runs(b):
    d = np.diff(np.r_[0, b.astype(int), 0]); return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def one(pid):
    vd = set(DAYS[DAYS.empatica_id == pid].date.astype(str))
    if len(vd) < 3:
        return None
    fm = os.path.join(P0, 'minute', f'{pid}.parquet')
    tz = pd.read_parquet(fm, columns=['tz']).tz.iloc[0]
    e5 = pd.read_parquet(os.path.join(P0, 'epoch5s', f'{pid}.parquet'), columns=['ts_utc', 'enmo', 'anglez'])
    f30 = pd.read_parquet(os.path.join(P0, 'features30s', f'{pid}.parquet'), columns=['ts_utc', 'worn', 'temp_c'])
    worn30 = f30.set_index('ts_utc').worn
    e5['worn'] = worn30.reindex(e5.ts_utc.dt.floor('30s')).fillna(False).values.astype(bool)
    e5['loc'] = e5.ts_utc.dt.tz_convert(tz).dt.tz_localize(None); e5['date'] = e5['loc'].dt.date.astype(str)
    e5 = e5[e5.date.isin(vd)]
    sp = NIGHTS[NIGHTS.empatica_id == pid]
    in_spt = np.zeros(len(e5), bool); t = e5.ts_utc.values
    for r in sp.itertuples():
        in_spt |= (t >= np.datetime64(r.spt_on.tz_convert(None))) & (t < np.datetime64(r.spt_off.tz_convert(None)))
    e5['in_spt'] = in_spt
    e5['dang'] = e5.anglez.diff().abs()
    D = []
    for dte, g in e5.groupby('date'):
        w = g[g.worn]
        if len(w) < 16 * 720:
            continue
        mg = w.enmo.values * 1000
        bins = np.arange(0, 4025, 25); cnt, _ = np.histogram(np.clip(mg, 0, 3999.9), bins); mins = cnt * 5 / 60; mid = bins[:-1] + 12.5
        k = mins > 0; ig = np.polyfit(np.log(mid[k]), np.log(mins[k]), 1)[0] if k.sum() >= 3 else np.nan
        srt = np.sort(mg)[::-1]
        mx = {f'm{x}_mg': srt[min(x * 12, len(srt)) - 1] for x in (60, 30, 15, 5)}
        m1 = w.set_index('loc').enmo.resample('1min').mean().dropna() * 1000
        act = (m1 >= 44.8).values; s, e = runs(act); s2, e2 = runs(~act)
        astp = 1 / np.mean(e - s) if len(s) else np.nan; satp = 1 / np.mean(e2 - s2) if len(s2) else np.nan
        # daytime sustained inactivity (van Hees 2015 SIB: >= 5 min with no 5-s arm-angle change > 5 deg),
        # outside the wrist sleep period, worn epochs; summed per day
        o = g[~g.in_spt & g.worn]
        s3, e3 = runs((o.dang <= 5).values)
        nap = sum(z - a for a, z in zip(s3, e3) if z - a >= 60) * 5 / 60          # 60 x 5-s epochs = 5 min
        nap30 = sum(z - a for a, z in zip(s3, e3) if z - a >= 360) * 5 / 60       # bouts >= 30 min only
        D.append(dict(date=dte, ig=ig, astp=astp, satp=satp, daytime_sib_min=nap, daytime_sib30_min=nap30, **mx))
    D = pd.DataFrame(D)
    if len(D) < 3:
        return None
    P = dict(empatica_id=pid, n_days_p2=len(D), **D.drop(columns='date').mean().to_dict())
    # steps
    st = pd.read_parquet(fm, columns=['ts_local', 'step_counts', 'wearing_detection_percentage'])
    st = st[st.ts_local.dt.date.astype(str).isin(vd) & (st.wearing_detection_percentage >= 50)]
    if len(st):
        g = st.groupby(st.ts_local.dt.date).step_counts
        P['steps_day'] = g.sum().mean(); P['peak30_cadence'] = g.apply(lambda s: s.nlargest(30).mean()).mean()
    # chronotype
    n = PN[PN.empatica_id == pid]
    if len(n) >= 3:
        fr = n.weekday.isin([4, 5])
        if fr.any() and (~fr).any():
            SOf, SDf, SDw = n[fr].onset_clock_h.mean(), n[fr].spt_h.mean(), n[~fr].spt_h.mean()
            SDweek = (5 * SDw + 2 * SDf) / 7
            P['msf_h'] = SOf + SDf / 2
            P['msfsc_h'] = P['msf_h'] if SDf <= SDw else SOf + SDweek / 2
    # temperature profile
    tc = f30[f30.worn].copy(); tc['loc'] = tc.ts_utc.dt.tz_convert(tz).dt.tz_localize(None)
    tc = tc[tc['loc'].dt.date.astype(str).isin(vd)].dropna(subset=['temp_c'])
    if len(tc) > 1000:
        prof = tc.groupby((tc['loc'].dt.hour * 6 + tc['loc'].dt.minute // 10)).temp_c.mean().reindex(range(144)).interpolate(limit_direction='both')
        circ = pd.concat([prof, prof]).reset_index(drop=True)
        r5 = circ.rolling(30).mean(); r10 = circ.rolling(60).mean()
        i5, i10 = int(r5.idxmax()), int(r10.idxmin())
        P.update(temp_m5=r5.max(), temp_l10=r10.min(), temp_m5_onset_h=((i5 - 29) % 144) / 6, temp_l10_onset_h=((i10 - 59) % 144) / 6)
        P['temp_ra'] = (P['temp_m5'] - P['temp_l10']) / (P['temp_m5'] + P['temp_l10'])
        th = (np.arange(144) + 0.5) / 6
        X = np.column_stack([np.ones(144), np.cos(2 * np.pi * th / 24), np.sin(2 * np.pi * th / 24)])
        b = np.linalg.lstsq(X, prof.values, rcond=None)[0]
        P['temp_cos_mesor'] = b[0]; P['temp_cos_acrophase_h'] = float((np.degrees(np.arctan2(b[2], b[1])) / 15) % 24)
    return P


if __name__ == '__main__':
    ids = sorted(DAYS.empatica_id.unique())
    with ProcessPoolExecutor(NW) as ex:
        P = pd.DataFrame([r for r in ex.map(one, ids) if r])
    P.to_csv(os.path.join(A4, 'phase2_blocks2to5_participant.csv'), index=False)
    pd.set_option('display.width', 200)
    print('participants', len(P))
    print(P.drop(columns=['empatica_id']).describe(percentiles=[.5]).T[['count', 'mean', '50%', 'min', 'max']].round(3).to_string())
