# Colorectal cancer scRNA-seq analysis with Scanpy

A reproducible single-cell RNA-seq analysis of colorectal tumors and adjacent non-malignant colon tissue using Scanpy.

## Objective

This project investigates cell-type composition and transcriptional-state differences between colorectal tumor and adjacent non-malignant colon tissue. It emphasizes reproducible data ingestion, transparent quality control, cell-type annotation, and patient-aware comparisons.

## Dataset

- Study: E-MTAB-8410
- Source: EMBL-EBI Single Cell Expression Atlas
- Assay: Single-cell RNA sequencing
- Tissue: Colorectal tumor, tumor border and adjacent non-malignant colon tissue
- Primary data format: Matrix Market sparse count matrix , Ensembl gene identifiers, and cell metadata

The original source files are not included in this repository. Download them from the EBI study directory and place them in `data/raw/`.

## Project structure

## Data ingestion

The EBI bundle includes:

- `E-MTAB-8410.aggregated_counts.mtx.gz`: unfiltered raw count matrix
- `E-MTAB-8410.aggregated_counts.mtx_rows.gz`: ordered Ensembl gene identifiers
- `E-MTAB-8410.aggregated_counts.mtx_cols.gz`: ordered cell identifiers
- `E-MTAB-8410.cell_metadata.tsv`: cell-level metadata

The source matrix has genes as rows and cells as columns. `src/data_ingest.py` transposes it to the AnnData convention: cells as rows and genes as columns.

The resulting AnnData object contains:

- `adata.X`: raw unfiltered counts
- `adata.layers["counts"]`: preserved raw-count copy
- `adata.obs`: cell metadata
- `adata.var`: Ensembl gene identifiers
- `adata.uns["data_provenance"]`: source and processing metadata

Open and run:

```bash
conda env create -f environment.yml
conda activate crc-scrnaseq

jupyter lab notebooks/01_build_raw_adata.ipynb
```

The ingestion workflow validates matrix dimensions, identifier uniqueness, metadata alignment, non-negative expression values, and integer-like counts.

## Data provenance and limitations

The raw count matrix is unfiltered. This project will define and document all downstream quality-control choices. Ensembl gene IDs are retained as the primary gene index; gene symbols will be added in a later annotation step.

## Quality control

Quality control was performed on the raw E-MTAB-8410 AnnData matrix before normalization, feature selection, dimensionality reduction, or clustering. Gene symbols were first mapped to Ensembl gene IDs using the GRCh37.87 GTF annotation so that mitochondrial genes could be identified from MT- gene symbols.
For each cell, the pipeline calculates four standard QC metrics:
- `n_genes_by_counts`: number of genes with nonzero expression
- `total_counts`: total expression/count value per cell
- `pct_counts_mt`: percentage of expression attributable to mitochondrial genes
- `pct_counts_ribo`: percentage of expression attributable to ribosomal genes
High mitochondrial expression combined with low library complexity can indicate damaged cells that have lost cytoplasmic RNA. Low detected-gene and count values can also represent empty droplets or low-quality cell barcodes. QC metrics were inspected globally and stratified by sample:`biosample_id`, because QC distributions differed substantially across libraries. High mitochondrial proportions are commonly used as a marker of low-quality cells, although suitable thresholds depend on tissue, protocol, and dataset context.

### Sample-level assessment
A per-sample QC summary was generated containing:
- Number of cells
- Median mitochondrial percentage
- Median detected genes
- Median total counts
- Fraction of cells with mitochondrial percentage above 20%
- Fraction of cells with fewer than 300 detected genes

Samples with median pct_counts_mt < 20% were classified as reference-good libraries. Libraries at or above this threshold were treated as potentially degraded and were not allowed to define their own QC baseline. This distinction was necessary because several libraries showed highly elevated mitochondrial distributions; applying within-sample outlier detection alone to those libraries could incorrectly treat broadly degraded cells as normal.

### Cell-level filtering
QC filtering uses two complementary approaches:
1. Reference-good samples: Cells were flagged using median absolute deviation (MAD) thresholds calculated independently within each good sample. Low total counts and low detected genes were evaluated on a log1p scale; high mitochondrial percentage was evaluated on its original percentage scale.
Absolute guardrails: A minimum detected-gene threshold was applied in addition to MAD-based filtering. The original study discarded droplets with fewer than 300 detected genes, so this threshold is retained as the primary empty-droplet filter for reproduction-oriented analysis. The study reports that droplets with number of detected genes below 300 were discarded.
2. Potential rescue from poor libraries: Cells from potentially degraded libraries were evaluated against thresholds derived from QC-passing cells from reference-good libraries. This prevents a sample-wide shift—such as a mitochondrial distribution centered at 70–80%—from making poor cells appear acceptable merely because they are not outliers within their own library.

The final QC mask retains:
- Cells from reference-good libraries that do not fail cell-level QC.
- Cells from potentially poor libraries only when they satisfy the reference-derived criteria and the global minimum-quality constraints.
All QC flags are retained in `adata.obs` to make filtering decisions auditable and to support per-library retention summaries.

### Doublets
The source study used DoubletFinder for doublet identification. DoubletFinder detects potential doublets by comparing observed cells with artificial doublets generated from combinations of cell-expression profiles.The supplied AnnData expression matrix contains substantive non-integer values. Therefore, the matrix cannot be safely treated as a raw UMI count matrix by rounding values for a raw-count-dependent doublet-calling workflow such as Scrublet.
For a faithful reproduction of the study’s doublet filtering, future work should use the authors’ original per-library count matrices and run DoubletFinder separately for each library, or recover author-provided doublet annotations if available. Until then, no independent doublet calls are made from the provided fractional-valued matrix; this limitation is recorded explicitly rather than applying a method outside its intended input assumptions.

## License and attribution

Data are provided by EMBL-EBI Expression Atlas under the source study’s terms. Cite E-MTAB-8410 and the associated publication when reusing the data.
https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-8410/
Joanito I, Wirapati P, Zhao N, Nawaz Z, Yeo G, Lee F, Eng CLP, Macalinao DC, Kahraman M, Srinivasan H, Lakshmanan V, Verbandt S, Tsantoulis P, Gunn N, Venkatesh PN, Poh ZW, Nahar R, Oh HLJ, Loo JM, Chia S, Cheow LF, Cheruba E, Wong MT, Kua L, Chua C, Nguyen A, Golovan J, Gan A, Lim WJ, Guo YA, Yap CK, Tay B, Hong Y, Chong DQ, Chok AY, Park WY, Han S, Chang MH, Seow-En I, Fu C, Mathew R, Toh EL, Hong LZ, Skanderup AJ, DasGupta R, Ong CJ, Lim KH, Tan EKW, Koo SL, Leow WQ, Tejpar S, Prabhakar S, Tan IB. Single-cell and bulk transcriptome sequencing identifies two epithelial tumor cell states and refines the consensus molecular classification of colorectal cancer. Nat Genet. 2022 Jul;54(7):963-975. doi: 10.1038/s41588-022-01100-4. Epub 2022 Jun 30. PMID: 35773407; PMCID: PMC9279158.
