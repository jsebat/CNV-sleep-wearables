"""Aim 1b report material: nights where the Withings mat adjudicates between headband and wrist.

Classes (from 13_triangulate, interval = headband end .. ORIGINAL wrist sleep-window end):
  'out of bed'                              mat: nobody on the mat -> wrist wrong (scored sleep while up)
  'in bed awake'                            mat: on the mat, awake -> wrist wrong (missed wake in bed)
  'headband stopped, still asleep in bed'   mat: asleep in bed -> headband truncated, wrist right
For each class the nights with the longest interval are shown (one per participant, N_EX per class),
each as a stacked panel: Dreem hypnogram | mat on-bed + sleep state | wrist P(sleep) original (12) and
retrained with mat out-of-bed examples (14) | wrist ENMO.  Vertical lines: headband end, mat out of bed,
original and retrained wrist window ends.  Also writes the class count table and the example list.
Usage: python -I 15_adjudication_examples.py <phase0_dir> <aim1_dir> [n_per_class]
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

P0, A1 = sys.argv[1:3]
NEX = int(sys.argv[3]) if len(sys.argv) > 3 else 3
OUT = os.path.join(A1, 'report', 'adjudication'); os.makedirs(OUT, exist_ok=True)
TR = pd.read_csv(os.path.join(A1, 'triangulation', 'triangulation_nightly.csv'))
for c in ['rec_start', 'rec_end', 'dreem_onset', 'dreem_final_wake', 'mat_in_bed', 'mat_out_bed', 'mat_final_sleep', 'wrist_on', 'wrist_off']:
    TR[c] = pd.to_datetime(TR[c], utc=True)
NEW = pd.read_csv(os.path.join(A1, 'retrain_mat', 'spt_nightly.csv'))
for c in ['spt_on', 'spt_off']:
    NEW[c] = pd.to_datetime(NEW[c], utc=True)
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
WT = pd.read_parquet(os.path.join(P0, 'withings_states.parquet'))
tzs = {}

# ---- class table ---------------------------------------------------------------------------------
M = TR[TR.has_mat == True]
tab = M.groupby('post_class').agg(nights=('recording_id', 'size'), participants=('empatica_id', 'nunique'),
                                  minutes=('post_min', 'sum'), median_min=('post_min', 'median')).reset_index()
verdict = {'out of bed': 'wrist wrong - scored sleep while out of bed', 'in bed awake': 'wrist wrong - missed wake in bed',
           'headband stopped, still asleep in bed': 'headband truncated - wrist right',
           'no post-headband interval': 'no disagreement to adjudicate'}
tab['verdict'] = tab.post_class.map(verdict)
tab.to_csv(os.path.join(OUT, 'adjudication_classes.csv'), index=False)
gen = M.pivot_table(index='genotype', columns='post_class', values='recording_id', aggfunc='count', fill_value=0)
gen.to_csv(os.path.join(OUT, 'adjudication_classes_by_genotype.csv'))


def new_window(r):
    c = NEW[(NEW.empatica_id == r.empatica_id) & (NEW.spt_on < r.dreem_final_wake + pd.Timedelta('6h')) & (NEW.spt_off > r.dreem_onset)]
    if not len(c):
        return pd.NaT, pd.NaT
    c = c.iloc[(np.minimum(c.spt_off, r.dreem_final_wake) - np.maximum(c.spt_on, r.dreem_onset)).values.argmax()]
    return c.spt_on, c.spt_off


ex = []
for cls in ['out of bed', 'in bed awake', 'headband stopped, still asleep in bed']:
    c = M[M.post_class == cls].sort_values('post_min', ascending=False).drop_duplicates('empatica_id').head(NEX)
    ex += [(cls, r) for r in c.itertuples()]
rows = []
STAGE_Y = {'W': 4, 'R': 3, 'N1': 2, 'N2': 1, 'N3': 0}
for k, (cls, r) in enumerate(ex, 1):
    pid = r.empatica_id
    if pid not in tzs:
        tzs[pid] = pd.read_parquet(os.path.join(P0, 'minute', f'{pid}.parquet'), columns=['tz']).tz.iloc[0]
    tz = tzs[pid]
    n_on, n_off = new_window(r)
    t0 = min(r.rec_start, r.mat_in_bed if pd.notna(r.mat_in_bed) else r.rec_start) - pd.Timedelta('1h')
    t1 = max(x for x in [r.rec_end, r.mat_out_bed, r.wrist_off, n_off] if pd.notna(x)) + pd.Timedelta('1h')
    old = pd.read_parquet(os.path.join(A1, 'deployable', f'scored_{pid}.parquet'))
    new = pd.read_parquet(os.path.join(A1, 'retrain_mat', f'scored_{pid}.parquet'))
    fe = pd.read_parquet(os.path.join(P0, 'features30s', f'{pid}.parquet'), columns=['ts_utc', 'enmo_mean'])
    sel = lambda df, col='ts_utc': df[(df[col] >= t0) & (df[col] <= t1)]
    d = sel(D[D.recording_id == r.recording_id]); old, new, fe = sel(old), sel(new), sel(fe)
    w = WT[(WT.study_id == r.study_id) & (WT.ts_utc >= t0) & (WT.ts_utc <= t1)]
    lt = lambda s: s.dt.tz_convert(tz).dt.tz_localize(None)
    lx = lambda x: x.tz_convert(tz).tz_localize(None) if pd.notna(x) else None
    fig, ax = plt.subplots(4, 1, figsize=(11, 7.2), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 1, 1.2, 0.8]))
    ax[0].step(lt(d.ts_utc), d.stage.map(STAGE_Y).astype(float), where='post', lw=1, color='#3b5b92')
    ax[0].set_yticks(range(5)); ax[0].set_yticklabels(['N3', 'N2', 'N1', 'REM', 'Wake']); ax[0].set_title('Dreem EEG headband', loc='left', fontsize=9)
    ax[1].scatter(lt(w.ts_utc), np.ones(len(w)) * 1.15, s=2, color='#999999', label='on mat')
    ax[1].step(lt(w.ts_utc), w.state.map({0: 0, 1: 0.6, 2: 0.3, 3: 0.9}), where='post', lw=1, color='#2e7d32')
    ax[1].set_yticks([0, 0.3, 0.6, 0.9, 1.15]); ax[1].set_yticklabels(['awake', 'deep', 'light', 'REM', 'on mat'], fontsize=7)
    ax[1].set_title('Withings mat (grey = someone on the mat)', loc='left', fontsize=9)
    ax[2].plot(lt(old.ts_utc), old.p_sleep.rolling(61, center=True, min_periods=20).median(), color='#c62828', lw=1.2, label='original wrist model')
    ax[2].plot(lt(new.ts_utc), new.p_sleep.rolling(61, center=True, min_periods=20).median(), color='#1565c0', lw=1.2, label='retrained with mat out-of-bed')
    ax[2].axhline(0.5, color='k', lw=0.5, ls=':'); ax[2].set_ylim(-0.05, 1.05); ax[2].set_ylabel('P(sleep)\n30-min median', fontsize=8)
    ax[2].legend(fontsize=7, loc='lower left'); ax[2].set_title('Wrist sleep probability', loc='left', fontsize=9)
    ax[3].plot(lt(fe.ts_utc), 1000 * fe.enmo_mean, lw=0.6, color='#555555'); ax[3].set_ylabel('ENMO mg', fontsize=8)
    ax[3].set_yscale('symlog', linthresh=5); ax[3].set_title('Wrist movement', loc='left', fontsize=9)
    for a in ax:
        for x, col, ls in ((r.rec_end, '#3b5b92', '-'), (r.mat_out_bed, '#2e7d32', '-'), (r.wrist_off, '#c62828', '--'), (n_off, '#1565c0', '--')):
            if pd.notna(x):
                a.axvline(lx(x), color=col, ls=ls, lw=1)
    fig.suptitle(f'Example {k}: {cls}  ->  {verdict[cls]}\n{r.study_id} ({r.genotype}), night of {lx(r.rec_start):%Y-%m-%d}; '
                 f'interval after headband {r.post_min:.0f} min, mat on-bed {100 * r.post_in_bed_frac:.0f} %',
                 fontsize=10, x=0.01, ha='left')
    fig.text(0.01, 0.005, 'vertical lines: blue solid = headband end | green = mat out of bed | red dashed = original wrist window end | '
             'blue dashed = retrained wrist window end', fontsize=7)
    plt.tight_layout(rect=(0, 0.02, 1, 0.94))
    fn = f'example_{k:02d}_{cls.split(",")[0].replace(" ", "_")}.png'
    fig.savefig(os.path.join(OUT, fn), dpi=160); plt.close(fig)
    rows.append(dict(example=k, file=fn, adjudication=cls, verdict=verdict[cls], study_id=r.study_id, genotype=r.genotype,
                     night=f'{lx(r.rec_start):%Y-%m-%d}', post_headband_min=r.post_min, mat_on_bed_frac=r.post_in_bed_frac,
                     mat_asleep_frac=r.post_asleep_frac,
                     orig_end_minus_mat_out=(r.wrist_off - r.mat_out_bed).total_seconds() / 60 if pd.notna(r.mat_out_bed) else np.nan,
                     new_end_minus_mat_out=(n_off - r.mat_out_bed).total_seconds() / 60 if pd.notna(n_off) and pd.notna(r.mat_out_bed) else np.nan))
E = pd.DataFrame(rows); E.to_csv(os.path.join(OUT, 'adjudication_examples.csv'), index=False)
pd.set_option('display.width', 220)
print(tab.to_string(index=False)); print(gen.to_string()); print(E.drop(columns=['file']).round(2).to_string(index=False))
