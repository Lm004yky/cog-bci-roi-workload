# ROI-Based EEG Spectral Markers of Cognitive Workload (COG-BCI reanalysis)

Analysis code and derived data for:
**Nurkhan & Gutoreva**, "ROI-Based EEG Spectral Markers of Cognitive Workload
and a Neural Recovery Index as a Proxy for Psychological Resilience."

## Data
Derived spectral features computed from the public COG-BCI dataset
(Hinss et al., 2023, Zenodo: https://zenodo.org/records/7413650). Raw EEG is
NOT redistributed here; only per-epoch/per-session derived metrics.

- `v2_all_epochs.csv` — per-epoch task features (N = 4,736 task epochs, 20 subjects)
- `v2_all_resting.csv` — per-epoch resting-state features
- `v2_nri_sessions.csv` — per-session Neural Recovery Index
- `v2_subject_workload_summary.csv` — per-subject summary

## Code
- `EEG_v2_cleanup.ipynb` — preprocessing + feature extraction + figures
- `code_fixes.py` — statistics (effect sizes, Holm/FDR correction, ICC + 95% CI,
  subject-level RM-ANOVA and cluster-robust sensitivity checks) and the Fig. 2 plot

## Reproducibility
- Python: <FILL: e.g. 3.11.x>   OS: <FILL: e.g. Ubuntu 22.04 / Colab>
- Package versions: see `requirements.txt`
- ICA: FastICA, 20 components, random_state = <FILL FROM YOUR PREPROCESSING CODE>,
  max_iter = <FILL>

To reproduce the statistics: `pip install -r requirements.txt`, then run the
BLOCK B/C/D cells of `code_fixes.py` on `v2_clean_epochs.csv` and
`v2_nri_sessions.csv`.

## License
<Add an MIT or CC-BY-4.0 license>
