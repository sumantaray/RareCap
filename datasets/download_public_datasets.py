#!/usr/bin/env python3
"""
Download public processed datasets used by RareCapsNet where a stable direct file URL
is available. Very large/raw datasets are intentionally not downloaded by default.
"""
from pathlib import Path
from urllib.request import urlretrieve

OUT = Path(__file__).resolve().parent / "downloaded"
OUT.mkdir(exist_ok=True)

FILES = {
    # PBMC68k filtered gene-barcode matrices (10x Genomics)
    "PBMC68k_filtered_gene_bc_matrices.tar.gz":
        "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz",

    # CBMC/CITE-seq processed RNA matrix (GEO GSE100866)
    "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",

    # Zeisel mouse brain processed expression matrix (GEO GSE60361)
    "GSE60361_C1-3005-Expression.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60361/suppl/GSE60361_C1-3005-Expression.txt.gz",

    # Darmanis processed per-cell files (GEO GSE67835)
    "GSE67835_RAW.tar":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67835/suppl/GSE67835_RAW.tar",
}

def download(name, url):
    dest = OUT / name
    if dest.exists():
        print(f"[skip] {dest.name} already exists")
        return
    print(f"[download] {name}")
    urlretrieve(url, dest)
    print(f"[done] {dest}")

def main():
    for name, url in FILES.items():
        try:
            download(name, url)
        except Exception as exc:
            print(f"[error] {name}: {exc}")

    print("\nAdditional datasets:")
    print("- Jurkat–293T: download the processed matrix from the 10x Genomics dataset page listed in README.md.")
    print("- Yan (GSE36552): GEO raw archive is extremely large; retrieve only the processed files actually used.")
    print("- Pollen (SRP041736): use NCBI SRA Toolkit (prefetch/fasterq-dump) or the processed matrix used by your analysis.")
    print("\nDo not commit multi-GB raw files to ordinary Git history.")

if __name__ == "__main__":
    main()
