"""Phase 0e: Withings Sleep Analyzer (under-mattress) exports -> sessions and per-minute states.

Inputs : /expanse/projects/sebat1/g2mh_data/withings/<site> SFARI/...  (folders and zips, possibly nested)
         Zips are extracted to <out>/withings_extracted/ (zip-slip guarded); the source is never modified.
         user.csv / account.csv (names, birthdates) are NOT read and NOT extracted.
Per export folder (any folder holding sleep.csv):
  sleep.csv              one row per bed session: from/to (in bed .. out of bed), light/deep/rem/awake s,
                         wake-ups, duration to sleep, duration to wake up. Headers are English or French.
  raw_bed_sleep-state.csv  start + lists of durations/values; values 0 awake, 1 light, 2 deep, 3 REM (README).
Study ID from the path: (MONT|TORO|PENN|UCLA)dddd[R2], or Montreal '211-<n>[_-]3ddd' -> MONT3ddd.
Exports overlap (Montreal per-day exports are cumulative) -> de-duplicated on (ID, session start) and
(ID, minute).
Outputs: <out>/withings_sessions.parquet, <out>/withings_states.parquet, <out>/withings_files.csv
Usage: python -I 05_withings.py <withings_root> <out_dir>
"""
import os, re, sys, json, zipfile, glob
import pandas as pd

SRC, OUT = sys.argv[1:3]
EX = os.path.join(OUT, 'withings_extracted'); os.makedirs(EX, exist_ok=True)
SKIP = {'user.csv', 'account.csv'}
IDRE = re.compile(r'(MONT|TORO|PENN|UCLA)(\d{4,5}(?:R2)?)')
MONTRE = re.compile(r'211-\d+[_-](3\d{3})')


def extract(z, dest):
    """extract a zip (and zips inside it) without writing outside dest and without identity files."""
    with zipfile.ZipFile(z) as zf:
        for m in zf.infolist():
            name = os.path.basename(m.filename)
            if m.is_dir() or name in SKIP:
                continue
            tgt = os.path.realpath(os.path.join(dest, m.filename))
            if not tgt.startswith(os.path.realpath(dest) + os.sep):
                continue
            os.makedirs(os.path.dirname(tgt), exist_ok=True)
            with zf.open(m) as r, open(tgt, 'wb') as w:
                w.write(r.read())
            if name.lower().endswith('.zip'):
                extract(tgt, os.path.splitext(tgt)[0]); os.remove(tgt)


for z in glob.glob(os.path.join(SRC, '**', '*.zip'), recursive=True):
    rel = os.path.relpath(z, SRC)
    dest = os.path.join(EX, os.path.splitext(rel)[0])
    if not os.path.isdir(dest):
        extract(z, dest)


def study_id(path):
    m = IDRE.search(path)
    if m:
        return m.group(1) + m.group(2)
    m = MONTRE.search(path)
    return 'MONT' + m.group(1) if m else None


COLS = {'from': 'start', 'début': 'start', 'to': 'end', 'fin': 'end', 'light (s)': 'light_s', 'léger (s)': 'light_s',
        'deep (s)': 'deep_s', 'profond (s)': 'deep_s', 'rem (s)': 'rem_s', 'awake (s)': 'awake_s', 'éveillé (s)': 'awake_s',
        'wake up': 'wakeups', 'réveils': 'wakeups', 'duration to sleep (s)': 'to_sleep_s', "durée d'endormissement (s)": 'to_sleep_s',
        'duration to wake up (s)': 'to_wake_s', 'durée de réveil (s)': 'to_wake_s', 'average heart rate': 'hr_mean',
        'fréquence cardiaque moyenne': 'hr_mean', 'night events': 'night_events'}
sess, states, files = [], [], []
for root in (SRC, EX):
    for f in glob.glob(os.path.join(root, '**', 'sleep.csv'), recursive=True):
        d = os.path.dirname(f); sid = study_id(os.path.relpath(d, root))
        files.append(dict(folder=os.path.relpath(d, root), source='extracted' if root == EX else 'original', study_id=sid))
        if not sid:
            continue
        s = pd.read_csv(f)
        s.columns = [COLS.get(c.strip().lower(), c) for c in s.columns]
        if 'start' in s and len(s):
            s['start_utc'] = pd.to_datetime(s.start, utc=True); s['end_utc'] = pd.to_datetime(s.end, utc=True)
            s['utc_offset_min'] = s.start.str[-6:].map(lambda o: (1 if o[0] == '+' else -1) * (int(o[1:3]) * 60 + int(o[4:6])))
            s['study_id'] = sid
            sess.append(s[[c for c in ['study_id', 'start_utc', 'end_utc', 'utc_offset_min', 'light_s', 'deep_s', 'rem_s', 'awake_s',
                                       'wakeups', 'to_sleep_s', 'to_wake_s', 'hr_mean', 'night_events'] if c in s]])
        fs = os.path.join(d, 'raw_bed_sleep-state.csv')
        if os.path.exists(fs):
            r = pd.read_csv(fs)
            rows = []
            for t0, du, va in zip(r.start, r.duration, r.value):
                t = pd.Timestamp(t0).tz_convert('UTC'); du = json.loads(du); va = json.loads(va)
                for a, b in zip(du, va):
                    for k in range(max(int(a) // 60, 1)):
                        rows.append((t + pd.Timedelta(minutes=k), b))
                    t += pd.Timedelta(seconds=a)
            if rows:
                x = pd.DataFrame(rows, columns=['ts_utc', 'state']); x['study_id'] = sid; states.append(x)
S = pd.concat(sess, ignore_index=True).drop_duplicates(['study_id', 'start_utc']).sort_values(['study_id', 'start_utc'])
ST = pd.concat(states, ignore_index=True).drop_duplicates(['study_id', 'ts_utc']).sort_values(['study_id', 'ts_utc'])
S.to_parquet(os.path.join(OUT, 'withings_sessions.parquet'), index=False)
ST.to_parquet(os.path.join(OUT, 'withings_states.parquet'), index=False)
F = pd.DataFrame(files); F.to_csv(os.path.join(OUT, 'withings_files.csv'), index=False)
print('export folders', len(F), '| without study ID', int(F.study_id.isna().sum()), F[F.study_id.isna()].folder.tolist()[:5])
print('participants', S.study_id.nunique(), '| sessions', len(S), '| state minutes', len(ST))
print('sessions per participant median', S.groupby('study_id').size().median())
print('UTC offsets used (min):', S.utc_offset_min.value_counts().to_dict())
S['hours'] = (S.end_utc - S.start_utc).dt.total_seconds() / 3600
print('session hours quantiles', S.hours.quantile([.05, .25, .5, .75, .95]).round(1).tolist())
