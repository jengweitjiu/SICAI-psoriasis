# SICAI: Stromal-Immune Coupled Attractor Index

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Computational framework for quantifying stromal-immune coupling architecture from spatial transcriptomics data, with application to psoriasis.

## Associated Publication

> **Geometric Stability Decomposition of Stromal-Immune Coupling Reveals Mast Cell Hub Architecture and Severity-Linked CD8⁺ TRM Attractor in Psoriasis**
>
> [Authors]. *Nature Communications* (2026). DOI: [to be added]

## Overview

SICAI integrates three metrics for each fibroblast–immune cell pair from Visium spatial transcriptomics:

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| **CS** (Coupling Strength) | Spatial Spearman ρ between fibroblast abundance and neighbourhood-averaged immune scores | Direction and magnitude of co-localisation |
| **CSp** (Coupling Specificity) | \|CS(i,j)\| / Σ\_k \|CS(i,k)\|, z-scored | Selectivity of partnership |
| **AS** (Attractor Stability) | 1 − CV across spatial windows | Spatial consistency |
| **SICAI** (Composite) | CS × CSp × AS | Biologically meaningful coupling rank |

**Geometric Stability Decomposition (GSD)** maps coupling states as attractors in a potential landscape, quantifying basin depth, coherence, and escape energy.

## Key Findings

- **Mast cell hub switch**: Healthy F4\_DP-centred coordination → lesional mast cell hub engaging all fibroblast subtypes
- **CD8⁺ TRM severity biomarker**: F2\_Universal↔CD8\_TRM coupling correlates with PASI (ρ = 0.73)
- **High-escape-energy attractor**: Lesional state has 2.8× escape energy vs healthy, explaining chronicity
- **CXCL12/CXCR4 therapeutic target**: 49% contribution to severity-linked TRM coupling

## Repository Structure

```
SICAI-psoriasis/
├── README.md
├── LICENSE                          # MIT License
├── CITATION.cff                     # Citation metadata for Zenodo
├── requirements.txt                 # Python dependencies
├── environment.yml                  # Conda environment
│
├── sicai/                           # Core SICAI module
│   ├── __init__.py
│   ├── coupling.py                  # CS, CSp, AS computation
│   ├── gsd.py                       # Geometric stability decomposition
│   ├── perturbation.py              # In silico L-R perturbation
│   └── utils.py                     # Helper functions
│
├── notebooks/                       # Analysis pipeline (Colab-ready)
│   ├── Part_A_Data_Download.py      # Download Visium + atlas data
│   ├── Part_B_Cell2location.py      # Fibroblast deconvolution
│   ├── Part_C_Immune_Scoring.py     # 13 immune population scoring
│   ├── Part_D_SICAI_Computation.py  # CS, CSp, AS, SICAI matrices
│   ├── Part_E_Clinical_Correlation.py  # PASI severity analysis
│   ├── Part_F_LR_Validation.py      # Ligand-receptor co-expression
│   ├── Part_G_GSD_Attractor.py      # Attractor landscape & basin metrics
│   ├── Part_H_Perturbation.py       # In silico perturbation analysis
│   ├── Part_I_Figure_Export.py      # Publication-quality figure generation
│   └── Part_J_Populate_Supp_Tables.py  # Fill supplementary tables
│
├── data/
│   └── immune_signatures.csv        # Curated gene signatures (Table S4)
│
└── figures/                         # Example outputs (not tracked in git)
    └── .gitkeep
```

## Installation

### Option 1: pip

```bash
git clone https://github.com/[username]/SICAI-psoriasis.git
cd SICAI-psoriasis
pip install -r requirements.txt
```

### Option 2: conda

```bash
git clone https://github.com/[username]/SICAI-psoriasis.git
cd SICAI-psoriasis
conda env create -f environment.yml
conda activate sicai
```

### Option 3: Google Colab (recommended for full pipeline)

Each notebook in `notebooks/` is designed to run on Colab with GPU (T4). Mount Google Drive and run sequentially:

```python
from google.colab import drive
drive.mount('/content/drive')
!pip install cell2location scanpy squidpy
```

## Quick Start

### Compute SICAI for your own data

```python
import scanpy as sc
import pandas as pd
from sicai.coupling import compute_cs, compute_csp, compute_as, compute_sicai

# Load your Visium AnnData with:
#   - adata.obsm['cell2location'] : fibroblast abundances (n_spots × n_subtypes)
#   - adata.obs[immune_cols]      : immune scores (n_spots × n_immune)
adata = sc.read_h5ad('your_visium.h5ad')

# Compute SICAI
cs_matrix = compute_cs(adata, fibroblast_key='cell2location', 
                        immune_keys=immune_cols, k_neighbours=6)
csp_matrix = compute_csp(cs_matrix)
as_matrix = compute_as(adata, fibroblast_key='cell2location',
                        immune_keys=immune_cols, window_size=500, n_windows=20)
sicai_matrix = compute_sicai(cs_matrix, csp_matrix, as_matrix)
```

### Geometric Stability Decomposition

```python
from sicai.gsd import fit_landscape, compute_basin_metrics

# cs_vectors: DataFrame (n_samples × n_pairs), rows = samples
landscape = fit_landscape(cs_vectors, bandwidth=0.4)
metrics = compute_basin_metrics(landscape, condition_labels)

print(f"Lesional escape energy: {metrics['lesional']['escape_energy']:.3f}")
print(f"Lesional coherence: {metrics['lesional']['coherence']:.3f}")
```

## Data Sources

| Dataset | Source | Accession |
|---------|--------|-----------|
| Psoriasis Visium | Castillo, Sidhu et al. (2023) | GEO: [GSE202011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202011) |
| Fibroblast atlas | Steele et al. (2023) | EBI: [S-BIAD2214](https://www.ebi.ac.uk/biostudies/studies/S-BIAD2214) |
| Processed SICAI data | This study | Zenodo: [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) |

## Processed Data (Zenodo)

The following processed datasets are deposited at Zenodo (DOI: 10.5281/zenodo.XXXXXXX):

- `CS_*.csv` — Coupling strength matrices (8 fibroblast × 13 immune) per condition
- `CSp_*.csv` — Coupling specificity matrices per condition
- `AS_*.csv` — Attractor stability matrices per condition
- `SICAI_*.csv` — Composite SICAI matrices per condition
- `sample_metrics.csv` — Per-sample coupling metrics with PASI
- `gsd_basin_metrics.csv` — GSD attractor basin metrics
- `transition_scores.csv` — Per-sample transition scores and PCA coordinates
- `perturbation_analysis.csv` — L-R perturbation contribution scores
- `Table_S1_all_pairs.csv` — Complete 104-pair SICAI results
- `Table_S2_LR_spatial.csv` — L-R spatial co-expression scores
- `Table_S3_sample_metrics.csv` — Per-sample clinical and computational metrics

## Software Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| scanpy | ≥1.9 | Single-cell/spatial analysis |
| cell2location | ≥0.1.3 | Spatial deconvolution |
| squidpy | ≥1.3 | Spatial neighbourhood analysis |
| scipy | ≥1.10 | Spearman correlation, KDE |
| scikit-learn | ≥1.2 | PCA, clustering |
| networkx | ≥3.0 | Coupling network visualisation |
| matplotlib | ≥3.7 | Figure generation |
| seaborn | ≥0.12 | Heatmaps |
| pandas | ≥1.5 | Data handling |
| numpy | ≥1.24 | Numerical computation |

## Pipeline Overview

```
Part A: Download Visium (GSE202011) + atlas (S-BIAD2214)
    ↓
Part B: Cell2location deconvolution → fibroblast abundances per spot
    ↓
Part C: Immune signature scoring → 13 populations per spot
    ↓
Part D: SICAI computation → CS, CSp, AS, SICAI matrices
    ↓
Part E: PASI correlation → severity biomarker identification
    ↓
Part F: L-R validation → spatial co-expression of 23 curated pairs
    ↓
Part G: GSD → attractor landscape, basin metrics, escape energy
    ↓
Part H: In silico perturbation → L-R contribution scores
    ↓
Part I: Figure export → 300 dpi TIFF/EPS/PDF panels
    ↓
Part J: Populate supplementary tables with pipeline results
```

## Reproducing Figures

Run Parts A–I sequentially on Colab (GPU recommended for Part B). Estimated runtime:

| Part | Time (T4 GPU) | Time (CPU) |
|------|---------------|------------|
| A | 10 min | 10 min |
| B | 2–4 hr | N/A (GPU required) |
| C | 15 min | 30 min |
| D | 20 min | 40 min |
| E | 5 min | 5 min |
| F | 10 min | 15 min |
| G | 10 min | 15 min |
| H | 15 min | 25 min |
| I | 5 min | 5 min |

## How to Cite

If you use SICAI in your research, please cite:

```bibtex
@article{SICAI2026,
  title={Geometric Stability Decomposition of Stromal-Immune Coupling Reveals 
         Mast Cell Hub Architecture and Severity-Linked CD8+ TRM Attractor 
         in Psoriasis},
  author={[Authors]},
  journal={Nature Communications},
  year={2026},
  doi={[to be added]}
}
```

And the software:

```bibtex
@software{SICAI_code2026,
  author={[Authors]},
  title={SICAI: Stromal-Immune Coupled Attractor Index},
  year={2026},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXXXX},
  url={https://github.com/[username]/SICAI-psoriasis}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

Correspondence: Tsen-Fang Tsai, M.D., Ph.D.  
Department of Dermatology, National Taiwan University Hospital  
Taipei, Taiwan
