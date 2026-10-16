# Multi-device wearable pipeline

Scripts are numbered in run order. Each script's docstring gives its inputs, outputs, definitions and usage.
All work on participant-level data happens on the cluster; only summary tables leave it.

## Phase 0: ingest and linkage

| Script | Step |
|---|---|
| `00_linkage.py` | Links wearable IDs to study GUIDs and CNV group (five groups: 16p del/dup, 22q del/dup, family control); flags repeat visits and test accounts |
| `01_minute_table.py` | Empatica per-minute vendor outputs (sleep detection, wear, pulse rate, etc.) into one table with local time zone |
| `02_raw_epochs.py` | Raw Avro → 5-s accelerometer epochs (sphere-fit calibration, ENMO, arm angle) and 30-s pulse, HRV, temperature and EDA |
| `03_dreem_epochs.py` | Headband hypnograms → UTC 30-s epochs |
| `04_daily_qc.py` | Wear time and valid days (≥ 16 h worn; ≥ 3 valid days per person) |
| `05_withings.py` | Mat sessions and per-minute in-bed / sleep states |
| `06_covariates.py` | Age at recording, BMI, medications classified by drug name (RxNorm lookup) |

## Aim 1: wrist band vs headband, triangulated with the mat

| Script | Step |
|---|---|
| `10_aim1_pairs.py`, `11_aim1_eval.py` | Paired wrist–headband nights; first-pass comparison of vendor, arm-angle and XGBoost sleep/wake |
| `12_aim1_deployable.py` | Deployable wrist sleep model (no headband-derived inputs) and 30-s feature table |
| `13_triangulate.py` | Three-device timing: clock alignment, lights-out, headband truncation |
| `14_aim1_retrain_mat.py` | Retraining with mat-certified out-of-bed epochs labelled wake |
| `15_adjudication_examples.py` | Worked nights where the mat decides between headband and wrist |
| `16_night_table.py`, `17_cn_models.R` | Night-level sleep table (TST, latency variants) and copy-number models |
| `18_figures_cn.R`, `19_figures_validation.py` | Figures |

## Aim 4: phenotype panel and copy-number effects

| Script | Step |
|---|---|
| `20_phenotypes.py` | Core wrist phenotypes per night, day and person |
| `21_aim4_models.R` | Mirror / U-shape / five-group mixed models; env `SEASON=1` (primary), `AGE_MIN/AGE_MAX` |
| `22_hr_checks.py`, `23_hr_models.R` | Heart-rate confounder checks (signal quality, still epochs, medication, BMI, heart defects) |
| `24_figures_aim4.R` | Figures (partial residuals from one reference model) |
| `25_phase2_hr_mat.py`, `26_phase2_models.R`, `27_phase2_more.py`, `28_eda_raw.py` | Second (exploratory) panel: heart-rate rhythm, mat physiology, activity structure, chronotype, temperature rhythm, EDA |
| `29_headband_core.py` | Headband per-recording metrics (depth index, REM latency, awakenings, latency from lights-out) |
| `30_master_table.py` | One row per participant across devices, plus a data-dictionary skeleton |

## Aim 2: XGBoost classifiers

| Script | Step |
|---|---|
| `35_aim2_xgb.py` | Nested family-grouped CV, in-fold confound removal, permutation null; env `NO_MEDS`, `SEASON`, `HB_SUBSET`, `FEATURE_SETS` (colon-separated) |
| `37_figures_aim2.py`, `38_aim2_features.py` | AUC figures; SHAP weights over family bootstraps joined to univariate effects |

## Aim 3: clinical features in carriers

| Script | Step |
|---|---|
| `40_aim3_outcomes.py` | Outcome table (diagnoses, FSIQ, Penn CNB accuracy and speed composites) |
| `41_aim3_assoc.R` | Association models with medication and one-carrier-per-family sensitivity |
| `42_aim3_predict.py` | Incremental prediction (wrist block permutation) |
| `44_aim3_stratified.R` | Sleep efficiency vs cognition by CNV group |

`run_*.sbatch` are the Slurm wrappers we used. Report and slide builders are kept internally because they embed
unpublished summary results.
