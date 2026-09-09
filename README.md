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

The original source files are not included in this repository. Download them from the EBI study directory (https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-8410/) and place them in `data/raw/`.

## Project Status

🚧 Work in progress

Completed:
- Data ingestion
- Quality-control preprocessing
- Gene-symbol attachment
- Normalization
- Clustering
- Cell-type annotation

Planned:
- Compare cell-type composition between tumor and adjacent non-malignant samples.
- Perform cell-type-specific pseudobulk differential expression analysis.
- Identify tumor-associated genes within malignant, immune, and stromal populations.
- Run pathway enrichment analysis on differentially expressed genes.
- Score tumor-related gene programs, including proliferation, hypoxia, and epithelial–mesenchymal transition.
- Characterize tumor-microenvironment changes across immune and stromal cell types.
- Explore cell–cell communication differences between tumor and adjacent tissue.
- Create publication-ready visualizations and summarize biological findings.


## Project structure

```text
├── data
│   ├── interim
│   │   ├── E-MTAB-8410_annotated.h5ad
│   │   ├── E-MTAB-8410_qc.h5ad
│   │   └── E-MTAB-8410_raw.h5ad
│   ├── raw
│   │   ├── E-MTAB-8410.aggregated_counts.mtx_cols.gz
│   │   ├── E-MTAB-8410.aggregated_counts.mtx_rows.gz
│   │   ├── E-MTAB-8410.aggregated_counts.mtx.gz
│   │   ├── E-MTAB-8410.cell_metadata.tsv
│   │   └── E-MTAB-8410.sdrf.txt
│   └── reference
│       └── Homo_sapiens.GRCh37.87.gtf.gz
├── environment.yml
├── figures
│   ├── qc_grouped_violin_metrics.png
│   ├── qc_total_counts_vs_genes.png
│   └── qc_violin_metrics.png
├── LICENSE
├── notebooks
│   ├── 01_create_raw_adata.ipynb
│   ├── 02_qc_preprocess.ipynb
│   ├── 03_normalization.ipynb
│   ├── 04_clustering.ipynb
│   └── 05_celltype_annotation.ipynb
├── README.md
├── results
└── src
    ├── attach_gene_symbols.py
    ├── ingest_data.py
    └── qc_preprocess.py
```


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

## Normalization
To preserve consistency with the original study, expression values were normalized using the same library-size normalization and log transformation strategy. Each cell was scaled to a total of 10,000 count-like expression units across all genes, followed by a natural-log transformation with a pseudocount of 1, log(1+x). This is equivalent to Seurat’s LogNormalize procedure with scale.factor = 10000 and was implemented in Scanpy using sc.pp.normalize_total(target_sum=10_000) followed by sc.pp.log1p(). The QC-retained input matrix, normalized counts, and log-normalized expression matrix were retained in AnnData layers for provenance and reproducibility.

## Clustering
Following principal component analysis, a k-nearest-neighbor graph was constructed using the first 30 principal components. A UMAP embedding was generated for visualization, and Leiden community detection was performed using the igraph implementation with a fixed random seed (random_state=0) to support reproducibility. Leiden clustering was evaluated at resolutions of 0.05, 0.1, and 0.5, allowing both broad cellular compartments and finer transcriptional substructure to be examined. Cluster resolution was selected based on UMAP structure, within-cluster consistency of automated cell-type predictions, and downstream marker-gene validation.

## Cell-type annotation
Cell types were initially assigned with CellTypist using the `Human_Colorectal_Cancer.pkl` reference model, which was selected to match the colorectal cancer tissue context of the dataset. CellTypist predictions were refined through majority voting and were retained as reference-based provisional labels together with classifier confidence scores. These predictions were compared with the study authors’ ontology labels, Leiden cluster composition, UMAP localization, and expression of canonical lineage markers. High-confidence, coherent populations included B cells, mast cells, enteric glial cells, stromal cells, myofibroblasts, pericytes, smooth muscle cells,and lymphatic  endothelial cells. Cells identified as immune cells were subsetted and re-annotated using `Immune_All_Low.pkl` immune reference model. Final annotations were based on concordance between the models and confidence scores. Where annotations were discordant or low confidence, broader immune lineages were retained.

## License and attribution

Data are provided by EMBL-EBI Expression Atlas under the source study’s terms. Cite E-MTAB-8410 and the associated publication when reusing the data.
https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-8410/

Lee HO, Hong Y, Etlioglu HE, Cho YB, Pomella V, Van den Bosch B, Vanhecke J, Verbandt S, Hong H, Min JW, Kim N, Eum HH, Qian J, Boeckx B, Lambrechts D, Tsantoulis P, De Hertogh G, Chung W, Lee T, An M, Shin HT, Joung JG, Jung MH, Ko G, Wirapati P, Kim SH, Kim HC, Yun SH, Tan IBH, Ranjan B, Lee WY, Kim TY, Choi JK, Kim YJ, Prabhakar S, Tejpar S, Park WY. Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer. Nat Genet. 2020 Jun;52(6):594-603. doi: 10.1038/s41588-020-0636-z. Epub 2020 May 25. PMID: 32451460.

Joanito I, Wirapati P, Zhao N, Nawaz Z, Yeo G, Lee F, Eng CLP, Macalinao DC, Kahraman M, Srinivasan H, Lakshmanan V, Verbandt S, Tsantoulis P, Gunn N, Venkatesh PN, Poh ZW, Nahar R, Oh HLJ, Loo JM, Chia S, Cheow LF, Cheruba E, Wong MT, Kua L, Chua C, Nguyen A, Golovan J, Gan A, Lim WJ, Guo YA, Yap CK, Tay B, Hong Y, Chong DQ, Chok AY, Park WY, Han S, Chang MH, Seow-En I, Fu C, Mathew R, Toh EL, Hong LZ, Skanderup AJ, DasGupta R, Ong CJ, Lim KH, Tan EKW, Koo SL, Leow WQ, Tejpar S, Prabhakar S, Tan IB. Single-cell and bulk transcriptome sequencing identifies two epithelial tumor cell states and refines the consensus molecular classification of colorectal cancer. Nat Genet. 2022 Jul;54(7):963-975. doi: 10.1038/s41588-022-01100-4. Epub 2022 Jun 30. PMID: 35773407; PMCID: PMC9279158.
