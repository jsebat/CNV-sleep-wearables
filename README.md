# CNV-sleep-wearables

At-home analysis of sleep, actigraphy and physiology from the joint analysis of multiple wearable devices,
applied to carriers of 16p11.2 and 22q11.2 deletions and duplications and their non-carrier family members
(G2MH / SFARI wearable study). Sebat Lab, UC San Diego.

Three devices record the same nights, each with a different blind spot, so they are used to check one another:

| Device | Signal | Used for |
|---|---|---|
| **Headband** (Dreem EEG) | 30-s sleep stages from EEG | reference sleep staging; sleep architecture (N3 depth, REM, latency) |
| **Wrist band** (Empatica EmbracePlus) | raw 64 Hz accelerometer, pulse (BVP + systolic peaks), skin temperature, electrodermal activity | sleep/wake on every night of ~2 weeks; activity, rest–activity rhythm, heart rate and HRV in sleep, EDA |
| **Mat** (Withings Sleep Analyzer) | under-mattress pressure, heart rate, breathing | exact in-bed / out-of-bed times; adjudicating headband–wrist disagreements; independent heart rate |

**This repository contains code only.** No participant-level data, phenotype tables, study identifiers or
results are included. Data are available to qualified researchers through the G2MH / SFARI data-access
procedures.

## Contents

| Folder | What it does |
|---|---|
| [`headband_sleep_architecture/`](headband_sleep_architecture/) | The headband-only analysis: Dreem hypnograms → QC → sleep metrics (TST, SOL, NREM depth index) → copy-number dose-response models (mirror and U-shape) → figures and report. R. |
| [`wearables/`](wearables/) | The multi-device analysis: raw wrist signal processing, a wrist sleep model trained on headband stages and corrected with the mat, a per-person phenotype panel, copy-number models, XGBoost classifiers and clinical associations. Python and R, run on a Slurm cluster. |
| [`scripts/id_scan.sh`](scripts/id_scan.sh) | Pre-commit check that no study identifier, GUID or credential has entered the repository. |

## Analysis aims (wearables)

1. **Validation.** How well does the wrist band recover sleep and sleep stages against the headband?
   An XGBoost epoch classifier is trained on headband stages and on mat-certified out-of-bed time, and
   evaluated on held-out families. The mat sets lights-out and adjudicates disagreements.
2. **Classifiers.** Can wearable phenotypes tell each CNV group from family members? Nested,
   family-grouped cross-validation, confounds regressed out in-fold, permutation tests, and SHAP feature weights
   compared against the univariate genotype effects.
3. **Clinical features.** Do wearable phenotypes relate to diagnoses and cognition (Penn CNB) in carriers,
   beyond CNV group? Association models with medication and relatedness checks, plus an incremental-prediction test.
4. **Phenotypes.** A standard per-person panel (sleep duration, timing and regularity; rest–activity rhythm;
   activity; heart rate, HRV and temperature in sleep), tested for mirror (linear in copy number) and U-shape
   effects with person and family random effects, adjusted for age, sex, site and season.

## Requirements

Python 3.11 (pandas, numpy, scipy, scikit-learn, xgboost, fastavro, pyarrow, matplotlib, joblib);
R 4.4 (lme4, lmerTest, ggplot2, sandwich, lmtest). Paths in the `run_*.sbatch` wrappers point at our cluster layout and
need to be edited for another system; the analysis scripts take their input and output locations as arguments.

## Citation

Manuscript in preparation. Please contact the Sebat Lab before reusing the code for publication.

## Licence

MIT (see [LICENSE](LICENSE)).
