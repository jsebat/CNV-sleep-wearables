"""Report figures + numbers for Aim 1 validation and triangulation.
Writes PNGs and summary_numbers.json to <out_dir>.
  val_tst_agreement.png   wrist TST - Dreem TST vs Dreem TST, original (12) vs retrained (14) model
  val_night_end.png       distribution of wrist night-end error vs the triangulated final awakening
  val_confusion4.png      4-class confusion (retrained model, out-of-fold, Dreem-labelled epochs), row %
Usage: python -I 19_figures_validation.py <phase0_dir> <aim1_dir> <out_dir>
"""
import os, sys, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

P0, A1, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#888888',
                     'axes.labelcolor': '#333333', 'xtick.color': '#555555', 'ytick.color': '#555555'})
ORIG, NEW, GREY = '#eb6834', '#2a78d6', '#8a8a86'
old = pd.read_csv(os.path.join(A1, 'deployable', 'spt_validation_nightly.csv')); old = old[(old.method == 'model') & (old.found == True)]
new = pd.read_csv(os.path.join(A1, 'retrain_mat', 'spt_validation_nightly.csv')); new = new[new.found == True]
tr = pd.read_csv(os.path.join(A1, 'triangulation', 'triangulation_nightly.csv'))
old = old.merge(tr[['recording_id', 'post_class']], on='recording_id', how='left')
# old end error vs triangulated reference
o2 = old.merge(new[['recording_id', 'end_err_triang', 'ref_extended']], on='recording_id', how='inner', suffixes=('', '_new'))
TRI = tr.set_index('recording_id')

fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
for a, d, est, ref, c, t in ((ax[0], old, 'est_tst', 'ref_tst', ORIG, 'Original wrist model'), (ax[1], new, 'est_tst', 'dreem_tst', NEW, 'Retrained with mat out-of-bed')):
    x = d[ref] / 60; yv = (d[est] - d[ref]) / 60
    a.scatter(x, yv, s=9, alpha=0.45, color=c, edgecolor='none')
    a.axhline(0, color='#333333', lw=0.8); m, s = yv.mean(), yv.std()
    for v, ls in ((m, '-'), (m - 1.96 * s, ':'), (m + 1.96 * s, ':')):
        a.axhline(v, color=c, lw=1.2, ls=ls)
    a.set_title(f'{t}\nmedian error {60 * yv.median():+.0f} min, bias {60 * m:+.0f} min', fontsize=10, loc='left')
    a.set_xlabel('Sleep measured by EEG headband (hours)')
ax[0].set_ylabel('Wrist minus headband (hours)'); ax[0].set_ylim(-6, 8)
plt.tight_layout(); fig.savefig(os.path.join(OUT, 'val_tst_agreement.png'), dpi=220); plt.close(fig)

fig, a = plt.subplots(figsize=(7.5, 3.4))
bins = np.arange(-300, 481, 20)
a.hist(o2.end_err_triang.clip(-300, 480), bins=bins, color=NEW, alpha=0.85, label='retrained', zorder=3)
a.hist(old.wake_err_min.clip(-300, 480), bins=bins, color=ORIG, alpha=0.55, label='original', zorder=2)
a.axvline(0, color='#333333', lw=0.8)
a.set_xlabel('Wrist night end minus true final awakening (minutes; clipped at -300 / +480)'); a.set_ylabel('Nights')
a.legend(frameon=False); plt.tight_layout(); fig.savefig(os.path.join(OUT, 'val_night_end.png'), dpi=220); plt.close(fig)

# confusion matrix from out-of-fold retrained predictions on Dreem epochs
D = pd.read_parquet(os.path.join(P0, 'dreem_epochs.parquet'))
D = D[(D.qc_pass_core == True) & ~D.stage.isin(['MT', 'U'])].assign(ts=lambda z: z.ts_utc.dt.round('30s'))
rows = []
for pid, d in D.groupby('dreem_id'):
    for cand in (pid, pid + 'R2'):
        f = os.path.join(A1, 'retrain_mat', f'scored_{cand}.parquet')
        if os.path.exists(f):
            s = pd.read_parquet(f, columns=['ts_utc', 'pred4']).dropna()
            rows.append(d.merge(s, left_on='ts', right_on='ts_utc', how='inner')[['stage', 'pred4']])
C = pd.concat(rows, ignore_index=True)
C['ref'] = C.stage.map({'W': 'Wake', 'N1': 'Light', 'N2': 'Light', 'N3': 'Deep', 'R': 'REM'})
C['est'] = C.pred4.map({0: 'Wake', 1: 'Light', 2: 'Deep', 3: 'REM'})
order = ['Wake', 'Light', 'Deep', 'REM']
M = pd.crosstab(C.ref, C.est).reindex(index=order, columns=order).fillna(0)
Mp = M.div(M.sum(1), axis=0) * 100
fig, a = plt.subplots(figsize=(4.8, 4.0))
a.imshow(Mp.values, cmap='Blues', vmin=0, vmax=100)
for i in range(4):
    for j in range(4):
        a.text(j, i, f'{Mp.values[i, j]:.0f}%', ha='center', va='center', color='white' if Mp.values[i, j] > 55 else '#222222', fontsize=10)
a.set_xticks(range(4), order); a.set_yticks(range(4), order)
a.set_xlabel('Wrist model says'); a.set_ylabel('EEG headband says')
a.spines[:].set_visible(False); plt.tight_layout(); fig.savefig(os.path.join(OUT, 'val_confusion4.png'), dpi=220); plt.close(fig)

# numbers for report/deck
r12 = json.load(open(os.path.join(A1, 'deployable', 'run_info.json'))); r14 = json.load(open(os.path.join(A1, 'retrain_mat', 'run_info.json')))
e11 = pd.read_csv(os.path.join(A1, 'results', 'epoch_binary.csv')).set_index('method')
S14 = pd.read_csv(os.path.join(A1, 'retrain_mat', 'spt_validation_summary.csv')).set_index('measure')
S12 = pd.read_csv(os.path.join(A1, 'deployable', 'spt_validation_summary.csv'))
adj = pd.read_csv(os.path.join(A1, 'report', 'adjudication', 'adjudication_classes.csv'))
num = dict(epoch_vendor_kappa=float(e11.loc['A_sleep', 'kappa']), epoch_angle_kappa=float(e11.loc['C_sleep', 'kappa']),
           epoch_vendor_spec=float(e11.loc['A_sleep', 'specificity']), epoch_angle_spec=float(e11.loc['C_sleep', 'specificity']),
           epoch_model_kappa=r14['epoch']['kappa'], epoch_model_sens=r14['epoch']['sensitivity'], epoch_model_spec=r14['epoch']['specificity'],
           epoch_model_acc=r14['epoch']['accuracy'], acc4=r14['epoch']['acc4'], kappa4=r14['epoch']['kappa4'], f1_4=r14['epoch']['F1_W_light_deep_REM'],
           out_of_bed_wake=r14['epoch']['out_of_bed_scored_wake'], nights_scored=r14['nights_scored'],
           onset_within30=float(S14.loc['sleep onset vs Dreem', 'within30']), onset_median=float(S14.loc['sleep onset vs Dreem', 'median']),
           end_within30_new=float(S14.loc['window end vs triangulated final sleep', 'within30']),
           end_median_new=float(S14.loc['window end vs triangulated final sleep', 'median']),
           tst_median_new=float(S14.loc['TST est - Dreem (all)', 'median']), tst_bias_new=float(S14.loc['TST est - Dreem (all)', 'bias']),
           tst_sd_new=float(S14.loc['TST est - Dreem (all)', 'sd']),
           end_median_old=float(S12[(S12.method == 'model') & (S12.measure == 'final awakening (all nights)')]['median'].iloc[0]),
           end_within30_old=float(S12[(S12.method == 'model') & (S12.measure == 'final awakening (all nights)')]['within30'].iloc[0]),
           tst_median_old=float(S12[(S12.method == 'model') & (S12.measure == 'TST (est - Dreem)')]['median'].iloc[0]),
           adjudication=adj[['post_class', 'nights', 'participants', 'minutes']].to_dict('records'),
           confusion_row_pct=Mp.round(1).to_dict())
json.dump(num, open(os.path.join(OUT, 'summary_numbers.json'), 'w'), indent=1)
print(json.dumps(num, indent=1)[:2500])
