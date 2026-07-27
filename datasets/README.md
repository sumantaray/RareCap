# RareCap datasets

This folder documents the public datasets used in the RareCapsNet study and provides reproducible download instructions.

## Important note

The original public datasets remain hosted by 10x Genomics, NCBI GEO, and NCBI SRA.
Large raw datasets should generally **not** be committed directly to a normal GitHub repository.
Instead, keep this folder in the repository with the manifest and download script below.

## Datasets

| Dataset | Public source | Accession / direct source |
|---|---|---|
| Jurkat–293T mixture | 10x Genomics | https://www.10xgenomics.com/datasets/50-percent-50-percent-jurkat-293-t-cell-mixture-1-standard-1-1-0 |
| PBMC68k | 10x Genomics | https://www.10xgenomics.com/datasets/fresh-68-k-pbm-cs-donor-a-1-standard-1-1-0 |
| CBMC/CITE-seq | NCBI GEO | GSE100866 |
| Zeisel mouse brain | NCBI GEO | GSE60361 |
| Yan | NCBI GEO | GSE36552 |
| Pollen | NCBI SRA | SRP041736 |
| Darmanis | NCBI GEO | GSE67835 |

## Suggested reproducibility wording

The public datasets analysed in this study are available from their original repositories.
This `datasets/` directory provides dataset identifiers, source URLs, and scripts/instructions for retrieving the data.

## Automated downloads

Run:

```bash
python download_public_datasets.py
```

This downloads the public processed files for which a stable direct URL is provided in the script.
Very large/raw archives are intentionally not downloaded by default.

For Pollen (SRP041736), install the NCBI SRA Toolkit and use `prefetch` / `fasterq-dump` as appropriate.

For the Jurkat–293T and PBMC68k 10x datasets, the 10x dataset pages provide the Cell Ranger output files.
