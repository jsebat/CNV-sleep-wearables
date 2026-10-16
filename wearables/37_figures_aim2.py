"""Aim 2 figures and summary tables (WRIST BAND, Empatica EmbracePlus).
Fig 1  AUC by comparison: A) feature set x age sample (main run); B) sensitivity runs (core, all ages).
Fig 2  TreeSHAP top features: 22q del vs control (core) and 16p del vs control (extended, wrist only).
Tables aim2_summary_main.csv, aim2_summary_sensitivity.csv
Usage: python -I 37_figures_aim2.py <aim2_dir> <out_dir>
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

A2, OUT = sys.argv[1:3]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
NAME = {('16pDel', 'Control'): '16p del vs control', ('16pDup', 'Control'): '16p dup vs control', ('22qDel', 'Control'): '22q del vs control',
        ('22qDup', 'Control'): '22q dup vs control', ('16pDel', '16pDup'): '16p del vs 16p dup', ('22qDel', '22qDup'): '22q del vs 22q dup'}
ORDER = list(NAME)
rd = lambda run: pd.read_csv(os.path.join(A2, run, 'aim2_performance.csv'))
P = pd.concat([rd('full'), rd('extended_wrist')]); P['run'] = 'main'
S = pd.concat([rd('full').assign(run='main'), rd('no_meds').assign(run='no medication'), rd('season').assign(run='+ season')])
keep = ['case', 'control', 'set', 'sample', 'n_case', 'n_ctrl', 'xgb_auc_mean', 'xgb_auc_lo', 'xgb_auc_hi', 'base_auc_mean', 'xgb_bal_acc_mean', 'xgb_cal_slope_mean', 'perm_p']
P[keep].round(3).to_csv(os.path.join(OUT, 'aim2_summary_main.csv'), index=False)
S[['run'] + keep].round(3).to_csv(os.path.join(OUT, 'aim2_summary_sensitivity.csv'), index=False)

def panel(ax, D, variants, title):
    y0 = np.arange(len(ORDER))[::-1]; off = np.linspace(-0.27, 0.27, len(variants))
    for (lab, sel, col, mk), o in zip(variants, off):
        for i, key in enumerate(ORDER):
            r = D[(D.case == key[0]) & (D.control == key[1])]
            r = r[np.logical_and.reduce([r[k] == v for k, v in sel.items()])]
            if not len(r): continue
            r = r.iloc[0]; y = y0[i] + o
            ax.plot([r.xgb_auc_lo, r.xgb_auc_hi], [y, y], color=col, lw=2, solid_capstyle='round')
            ax.plot(r.xgb_auc_mean, y, mk, color=col, ms=6.5, mec='white', mew=0.8, label=lab if i == 0 else None)
            ax.plot(r.base_auc_mean, y, '|', color='#777777', ms=7, mew=1.4, label='age/sex/site only' if (i == 0 and o == off[0]) else None)
            if r.perm_p < 0.05: ax.text(max(r.xgb_auc_hi, r.base_auc_mean) + 0.012, y, f'p={r.perm_p:.3f}', va='center', fontsize=7.2, color=col)
    ax.axvline(0.5, color='#999999', lw=0.8, ls='--'); ax.set_xlim(0.2, 1.02)
    ax.set_yticks(np.arange(len(ORDER))[::-1]); ax.set_yticklabels([NAME[k] for k in ORDER])
    ax.set_xlabel('AUC (point = mean of 20 repeats; bar = 2.5–97.5 % of repeats)'); ax.set_title(title, loc='left', fontsize=10, fontweight='bold')
    ax.legend(loc='lower left', fontsize=7.3, frameon=False)

fig, axs = plt.subplots(1, 2, figsize=(12.2, 5.0), sharey=True)
panel(axs[0], P, [('core, all ages', dict(set='core', sample='all'), '#1f4e79', 'o'),
                  ('core, age 7–25', dict(set='core', sample='age7to25'), '#6a9bc3', 's'),
                  ('extended (wrist only), all ages', dict(set='extended_wrist', sample='all'), '#c0504d', 'D')],
      'A  Feature set and age range (main run)')
panel(axs[1], S, [('main', dict(run='main', set='core', sample='all'), '#1f4e79', 'o'),
                  ('no medication', dict(run='no medication', set='core', sample='all'), '#9b6a00', 's'),
                  ('+ season', dict(run='+ season', set='core', sample='all'), '#4f7f3a', 'D')],
      'B  Sensitivity runs (core set, all ages)')
fig.suptitle('Wrist band (Empatica EmbracePlus): how well participant-level phenotypes separate each CNV group', x=0.01, ha='left', fontsize=10.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'aim2_auc.png'), dpi=200); plt.close(fig)

LBL = {'wr_rmssd_sleep': 'HRV in sleep (RMSSD)', 'wr_hr_sleep_mean': 'mean HR in sleep', 'wr_hr_sleep_min': 'lowest HR in sleep', 'wr_IV': 'intradaily variability',
       'wr_tst_h': 'total sleep time', 'wr_spt_h': 'sleep period', 'wr_mvpa_min': 'moderate-vigorous activity', 'wr_temp_sleep_mean': 'wrist temp in sleep',
       'wr_m10_mg': 'most-active 10 h level', 'wr_p2_eda_peaks_per_h_day': 'EDA peaks per hour, day', 'wr_p2_m15_mg': 'most-active 15 min level',
       'wr_p2_peak30_cadence': 'peak 30-min cadence', 'wr_p2_hr_day_still': 'daytime HR when still', 'wr_p2_m60_mg': 'most-active 60 min level',
       'wr_p2_astp': 'active-to-sedentary transitions', 'wr_social_jetlag_h': 'social jet lag', 'wr_tst_sd_h': 'SD of sleep time', 'wr_SRI': 'Sleep Regularity Index',
       'wr_p2_m30_mg': 'most-active 30 min level', 'wr_m10_onset_h': 'most-active 10 h onset', 'wr_light_h': 'light activity', 'wr_l5_mg': 'least-active 5 h level',
       'wr_sed_h': 'sedentary time', 'wr_cos_mesor': 'cosinor mesor', 'wr_temp_amp': 'wrist temp range', 'wr_IS': 'interdaily stability',
       'wr_cos_acrophase_h': 'activity peak time', 'wr_p2_eda_scl_day': 'skin conductance level, day', 'wr_p2_daytime_sib30_min': 'daytime inactivity ≥30 min',
       'wr_p2_temp_cos_acrophase_h': 'temperature rhythm peak time', 'wr_p2_ig': 'intensity gradient', 'wr_p2_hr_dip_pct': 'night-time HR dip',
       'wr_p2_temp_l10_onset_h': 'temperature L10 onset', 'wr_p2_eda_scl_sleep': 'skin conductance level, sleep', 'wr_p2_hr_cos_amp': 'HR rhythm amplitude',
       'wr_p2_hr_cos_mesor': 'HR rhythm mesor', 'wr_onset_clock_h': 'sleep onset time', 'wr_p2_temp_cos_mesor': 'temperature rhythm mesor',
       'wr_RA': 'relative amplitude', 'wr_midpoint_sd_h': 'SD of sleep midpoint'}
def shp(run, case, ctrl, st):
    # bootstrap-median weights (38_aim2_features.py)
    E = pd.read_csv(os.path.join(A2, 'features', 'aim2_feature_weights.csv')).rename(columns={'importance': 'mean_abs_shap'})
    return E[(E.case == case) & (E.control == ctrl) & (E.set == st)].sort_values('rank').head(8)
fig, axs = plt.subplots(1, 2, figsize=(12.2, 3.9))
for ax, (run, case, ctrl, st, ttl) in zip(axs, [('full', '22qDel', 'Control', 'core', 'A  22q del vs control (core set)'),
                                               ('extended_wrist', '16pDel', 'Control', 'extended_wrist', 'B  16p del vs control (extended, wrist only)')]):
    t = shp(run, case, ctrl, st).iloc[::-1]
    col = np.where(t.shap_direction > 0, '#c0504d', '#2b6ca3')
    ax.barh([LBL.get(f, f) for f in t.feature], t.mean_abs_shap, color=col, height=0.62)
    for yy, (v, ok) in enumerate(zip(t.mean_abs_shap, t.sign_agrees)):
        ax.text(v + 0.006, yy, '' if ok != 0 else 'sign differs from univariate effect', va='center', fontsize=7, color='#666666')
    ax.set_xlabel('weight: median over 100 bootstraps of mean |SHAP| (log-odds)'); ax.set_title(ttl, loc='left', fontsize=10, fontweight='bold')
from matplotlib.patches import Patch
axs[0].legend(handles=[Patch(color='#c0504d', label='higher value → more carrier-like'), Patch(color='#2b6ca3', label='lower value → more carrier-like')],
              loc='lower right', fontsize=7.5, frameon=False)
fig.suptitle('Wrist band (Empatica EmbracePlus): features driving each model (TreeSHAP, 100 family bootstrap refits)', x=0.01, ha='left', fontsize=10.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'aim2_shap.png'), dpi=200); plt.close(fig)
# Fig 3: feature-weight map across comparisons (extended, wrist only; bootstrap-median weights from 38_aim2_features.py)
FW = os.path.join(A2, 'features', 'aim2_feature_weights.csv')
if os.path.exists(FW):
    W = pd.read_csv(FW); W = W[W.set == 'extended_wrist']
    W['imp_rel'] = W.importance / W.groupby(['case', 'control']).importance.transform('max')
    feats = list(dict.fromkeys(W[W['rank'] <= 6].sort_values(['case', 'control', 'rank']).feature))
    cols = ORDER
    Z = np.full((len(feats), len(cols)), np.nan); Pv = np.full_like(Z, np.nan)
    for j, (c, k) in enumerate(cols):
        t = W[(W.case == c) & (W.control == k)].set_index('feature')
        for i, f in enumerate(feats):
            if f in t.index: Z[i, j] = t.loc[f, 'imp_rel'] * t.loc[f, 'shap_direction']; Pv[i, j] = t.loc[f, 'univ_p']
    fig, ax = plt.subplots(figsize=(8.6, 0.36 * len(feats) + 1.9))
    im = ax.imshow(Z, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(len(feats)):
        for j in range(len(cols)):
            if np.isfinite(Pv[i, j]) and Pv[i, j] < 0.05: ax.text(j, i, '●', ha='center', va='center', fontsize=6.5, color='black')
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([LBL.get(f, f.replace('wr_p2_', '').replace('wr_', '')) for f in feats], fontsize=8)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([NAME[k] for k in cols], rotation=30, ha='right', fontsize=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.6); cb.set_label('weight relative to the top feature\n(red: higher value → more like the first group)', fontsize=7.5)
    ax.set_title('Wrist band: feature weights by comparison (extended set, wrist only)\n● = univariate genotype/dosage effect p < 0.05 (Aim 4 / second-panel models)',
                 loc='left', fontsize=9.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'aim2_feature_map.png'), dpi=200); plt.close(fig)
print('figures written to', OUT)
