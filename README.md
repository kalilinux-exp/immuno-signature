[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21074547.svg)](https://doi.org/10.5281/zenodo.21074547)
# ImmunoScore

A CD8 T-cell gene signature that scores anti-PD-1 (checkpoint immunotherapy)
response from single-cell tumor data. The signature was built on melanoma and
tested on basal cell carcinoma, a cancer it was never trained on.

**Result:** patient-level cross-cancer AUC **0.77** (n = 11 patients,
permutation p = 0.08, 95% CI [0.40, 1.00]). The signal is suggestive and
hypothesis-generating, not statistically significant at this sample size. The
robust finding is the recovered biology: a stem-like / memory CD8 program
(TCF7, CCR7, IL7R, CD28) already linked to checkpoint response.

Research and education use only. Not a medical device.

## What's here

| Path | Purpose |
|------|---------|
| `01_learn_pipeline.py` | Learn and verify the single-cell workflow on a small built-in dataset |
| `02_real_data.py` | Build the responder signature from melanoma (GSE120575) |
| `03_validate.py` | Score the signature on basal cell carcinoma (GSE123813); cell + patient AUC |
| `04_pseudobulk.py` | Patient-level (pseudobulk) signature, the more rigorous version |
| `05_figures.py` | Volcano plot and CD8 UMAP |
| `06_rigor.py` | Bootstrap confidence interval and permutation test |
| `scpipe.py` | Shared, unit-tested helper functions |
| `test_pipeline.py` | Unit tests (run with plain Python, no framework) |
| `app.py` | ImmunoScore Streamlit tool |
| `docs/` | Project website (GitHub Pages) |
| `*_signature.txt` | The saved gene signatures |
| `figures/` | Generated plots |

## Setup

Python 3.11 or newer (developed on 3.14).

```bash
py -3 -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

## Data (not included — download it)

The single-cell datasets are large, so they are not committed. Download both
from [TISCH2](http://tisch.comp-genomics.org) and place them in `data/`:

1. **Melanoma** — search `SKCM_GSE120575_aPD1aCTLA4`. From its Download menu, get
   *Single-cell expression matrix* (`.h5`) and *Meta information* (`.tsv`).
2. **Basal cell carcinoma** — browse Skin, BCC, dataset `GSE123813`. Download the
   same two files.

Your `data/` folder should contain:
```
data/
  SKCM_GSE120575_aPD1aCTLA4_expression.h5
  SKCM_GSE120575_aPD1aCTLA4_CellMetainfo_table.tsv
  BCC_GSE123813_aPD1_expression.h5
  BCC_GSE123813_aPD1_CellMetainfo_table.tsv
```

## Reproduce

```bash
.venv\Scripts\python 02_real_data.py     # melanoma signature -> response_signature.txt
.venv\Scripts\python 03_validate.py      # cross-cancer validation
.venv\Scripts\python 04_pseudobulk.py    # patient-level signature, AUC 0.77
.venv\Scripts\python 06_rigor.py         # confidence interval + permutation test
.venv\Scripts\python test_pipeline.py    # unit tests
.venv\Scripts\streamlit run app.py       # interactive tool
```

## Limitations

- **Small sample.** Eleven validation patients. Wide confidence interval, p = 0.08.
- **Public reanalysis.** Built from published datasets, not new patient data.
- **Possible confounds.** Two datasets from different labs can carry batch effects;
  the signature may partly reflect general T-cell health rather than response alone.

## Data sources

- Sade-Feldman et al., *Cell*, 2018 (melanoma, GSE120575)
- Yost et al., *Nature Medicine*, 2019 (basal cell carcinoma, GSE123813)
- Preprocessed via TISCH2 (Sun et al., *Nucleic Acids Research*, 2021)

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 Kalixte Petrof.
