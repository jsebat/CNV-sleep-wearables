"""Phase 0d: Dreem hypnograms -> one 30-s epoch table on the UTC time grid.

Filename: SFARI_Dreem_Prod_Hypno_<ID>_<YYYY-MM-DD>T<HH_MM_SS><+/-HH_MM>__<n>.txt
The filename gives the recording start in local time with its UTC offset; epochs are
consecutive 30-s rows. Stage labels: SLEEP-S0 (wake), S1, S2, S3, SLEEP-REM, SLEEP-MT
(movement), anything else -> unscored.
Night QC flags are joined from the sleep-architecture per-recording table, so the same
782/761 nights can be selected.

Usage: python -I 03_dreem_epochs.py <hypnogram_dir> <per_recording_metrics_with_qc.csv> <out_dir>
"""
import os, re, sys, glob
import pandas as pd

HD, QC, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
STAGE = {'SLEEP-S0': 'W', 'SLEEP-S1': 'N1', 'SLEEP-S2': 'N2', 'SLEEP-S3': 'N3', 'SLEEP-REM': 'R', 'SLEEP-MT': 'MT'}
pat = re.compile(r'Hypno_(?P<id>[A-Z]+\d+)\s*_(?P<d>\d{4}-\d{2}-\d{2})T(?P<t>\d{2}_\d{2}_\d{2})(?P<off>[+-]\d{2}_\d{2})__(?P<n>\d+)')
rows, bad = [], []
for f in sorted(glob.glob(os.path.join(HD, '*.txt'))):
    b = os.path.basename(f)
    m = pat.search(b)
    if not m:
        bad.append(b); continue
    start = pd.Timestamp(f"{m['d']}T{m['t'].replace('_', ':')}{m['off'].replace('_', ':')}")
    lines = open(f, encoding='utf-8', errors='replace').read().splitlines()
    try:
        h = next(i for i, l in enumerate(lines) if l.startswith('Sleep Stage'))
    except StopIteration:
        bad.append(b); continue
    st = [l.split('\t')[0].strip() for l in lines[h + 1:] if l.strip()]
    rid = os.path.splitext(b)[0]
    rows.append(pd.DataFrame({'dreem_id': m['id'], 'recording_id': rid, 'epoch': range(len(st)),
                              'ts_utc': start.tz_convert('UTC') + pd.to_timedelta([30 * i for i in range(len(st))], unit='s'),
                              'utc_offset_min': int(start.utcoffset().total_seconds() / 60),
                              'stage': [STAGE.get(s, 'U') for s in st]}))
E = pd.concat(rows, ignore_index=True)
q = pd.read_csv(QC)
key = 'recording_id' if 'recording_id' in q.columns else q.columns[0]
qc_cols = [c for c in q.columns if c.startswith('excl_') or c == 'qc_pass_core']
E = E.merge(q[[key] + qc_cols].drop_duplicates(key), left_on='recording_id', right_on=key, how='left')
E.to_parquet(os.path.join(OUT, 'dreem_epochs.parquet'), index=False)
n = E.groupby('recording_id').agg(id=('dreem_id', 'first'), qc=('qc_pass_core', 'first'))
print('recordings', len(n), 'participants', n.id.nunique(), 'qc-pass', int(n.qc.fillna(False).astype(bool).sum()),
      'unmatched-to-QC', int(n.qc.isna().sum()), 'unparsed files', len(bad))
print(E.stage.value_counts().to_string())
