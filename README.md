# RareCapsNet

RareCapsNet is an interpretable capsule-network framework for identifying rare
cell populations in single-cell transcriptomic data and examining
capsule-associated gene signals.

This repository currently contains a partial reproducibility package for the
available RareCapsNet analyses.

## Repository contents

- `notebooks/`: Jurkat and simulation analysis notebooks. Stored outputs have
  been removed to keep the repository compact.
- `scripts/`: Python and R analysis scripts.
- `datasets/`: Public-dataset manifest and download instructions.
- `data/simulated/`: Included S1 and S3 simulated expression matrices and
  labels.
- `data/processed/`: Small processed data used by the analysis notebooks.
- `results/jurkat/`: Jurkat coupling and selected-gene results.
- `figures/`: Workflow and simulation figures in PDF and SVG formats.

## Python environment

```bash
conda create -n rarecap python=3.8
conda activate rarecap
pip install -r requirements.txt
```

The notebooks use PyTorch, NumPy, pandas, scikit-learn, and Matplotlib. The R
analysis uses packages including `ComplexHeatmap`, `RColorBrewer`, and
`ggplot2`.

## Data

Public datasets remain hosted by their original repositories. See
[`datasets/README.md`](datasets/README.md) for identifiers and download
instructions. The included S1 and S3 simulation files are located under
`data/simulated/`.

Some notebooks retain the original research filenames or input-path
assumptions. Adjust the input paths to the corresponding local dataset before
running an analysis.
