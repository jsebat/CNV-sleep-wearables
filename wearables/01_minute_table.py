"""Phase 0b: assemble Empatica per-minute digital biomarkers into one table per participant.

Input  : <root>/1/1/participant_data/<UTC date>/<ID>-<serial>/digital_biomarkers/aggregated_per_minute/*.csv
         <root>/1/1/metadata/1-1-<ID>_metadata.csv  (timezone_location, time_offset)
Output : <out>/minute/<ID>.parquet - one row per minute (UTC), all biomarkers wide,
         plus local time (IANA zone from metadata; DST-correct) and the vendor
         missing-value reasons for wear and sleep detection.
Usage  : python -I 01_minute_table.py <root> <out_dir> [n_workers]
"""
import os, sys, glob, re
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

ROOT, OUT = sys.argv[1], sys.argv[2]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
PD = os.path.join(ROOT, '1', '1', 'participant_data')
MD = os.path.join(ROOT, '1', '1', 'metadata')
os.makedirs(os.path.join(OUT, 'minute'), exist_ok=True)

SITE_TZ = {'MONT': 'America/Toronto', 'TORO': 'America/Toronto', 'PENN': 'America/New_York', 'UCLA': 'America/Los_Angeles'}
ALIAS = {'US/Pacific': 'America/Los_Angeles', 'US/Central': 'America/Chicago', 'US/Mountain': 'America/Denver',
         'US/Eastern': 'America/New_York', 'America/Montreal': 'America/Toronto'}
# value column(s) kept from each biomarker file
KEEP = {'accelerometers-std': ['accelerometers_std_g'],
        'actigraphy-counts': ['counts_x_axis', 'counts_y_axis', 'counts_z_axis', 'vector_magnitude'],
        'activity-classification': ['activity_class'], 'activity-counts': ['activity_counts'],
        'activity-intensity': ['activity_intensity'], 'body-position': ['body_position_left', 'body_position_right'],
        'eda': ['eda_scl_usiemens'], 'met': ['met'], 'prv': ['prv_rmssd_ms'], 'pulse-rate': ['pulse_rate_bpm'],
        'respiratory-rate': ['respiratory_rate_brpm'], 'sleep-detection': ['sleep_detection_stage'],
        'step-counts': ['step_counts'], 'temperature': ['temperature_celsius'],
        'wearing-detection': ['wearing_detection_percentage']}
REASON = {'wearing-detection': 'wear_missing_reason', 'sleep-detection': 'sleep_missing_reason'}


def participant_tz(pid):
    base = re.sub(r'R2$', '', pid)
    tz = None
    for f in (os.path.join(MD, f'1-1-{pid}_metadata.csv'), os.path.join(MD, f'1-1-{base}_metadata.csv')):
        if os.path.exists(f):
            m = pd.read_csv(f, dtype=str)
            z = m.get('timezone_location', pd.Series(dtype=str)).dropna()
            if len(z):
                tz = z.value_counts().index[0]
                break
    src = 'metadata' if tz else 'site default'
    tz = ALIAS.get(tz, tz) or SITE_TZ.get(re.match(r'[A-Z]+?(?=\d|TEST)', pid).group(0), 'UTC')
    return tz, src


def build(pid):
    dirs = sorted(glob.glob(os.path.join(PD, '*', f'{pid}-*', 'digital_biomarkers', 'aggregated_per_minute')))
    frames = []
    for d in dirs:
        parts = None
        for kind, cols in KEEP.items():
            fs = glob.glob(os.path.join(d, f'*_{kind}.csv'))
            if not fs:
                continue
            x = pd.read_csv(fs[0], usecols=lambda c: c in ['timestamp_unix', 'missing_value_reason'] + cols)
            if kind in REASON:
                x = x.rename(columns={'missing_value_reason': REASON[kind]})
            else:
                x = x.drop(columns=['missing_value_reason'], errors='ignore')
            parts = x if parts is None else parts.merge(x, on='timestamp_unix', how='outer')
        if parts is not None:
            parts['device_sn'] = os.path.basename(os.path.dirname(os.path.dirname(d))).split('-', 1)[1]
            frames.append(parts)
    if not frames:
        return pid, 0, None, None
    df = pd.concat(frames).drop_duplicates('timestamp_unix').sort_values('timestamp_unix')
    tz, src = participant_tz(pid)
    df['ts_utc'] = pd.to_datetime(df.timestamp_unix, unit='ms', utc=True)
    loc = df.ts_utc.dt.tz_convert(tz)
    df['ts_local'] = loc.dt.tz_localize(None)
    df['utc_offset_min'] = (loc.apply(lambda t: t.utcoffset().total_seconds() / 60)).astype('int16')
    df['empatica_id'] = pid
    df['tz'] = tz
    df.to_parquet(os.path.join(OUT, 'minute', f'{pid}.parquet'), index=False)
    return pid, len(df), tz, src


if __name__ == '__main__':
    ids = sorted({os.path.basename(p).split('-', 1)[0] for p in glob.glob(os.path.join(PD, '*', '*'))})
    with ProcessPoolExecutor(NW) as ex:
        res = list(ex.map(build, ids))
    s = pd.DataFrame(res, columns=['empatica_id', 'n_minutes', 'tz', 'tz_source'])
    s.to_csv(os.path.join(OUT, 'minute_build_summary.csv'), index=False)
    print(s.tz_source.value_counts().to_string())
    print('participants', len(s), 'minutes', s.n_minutes.sum())
