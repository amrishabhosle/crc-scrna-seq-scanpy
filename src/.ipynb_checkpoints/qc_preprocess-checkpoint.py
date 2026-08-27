from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt


def calculate_qc_metrics(
    adata: ad.AnnData
) -> None:
    """Calculate total counts, total genes, and mitochondrial genes."""

    if ("gene_symbol" not in adata.var.columns and adata.var_names.str.startswith("ENSG").all()):
        raise ValueError("Gene symbols not found. Run attach_gene_symbol() first.")

    if not adata.var_names.is_unique:
        raise ValueError("Gene identifiers in adata.var_names must be unique.")

    gene_symbols = adata.var["gene_symbol"].astype(str)

    # Mitochondrial genes
    adata.var["mt"] = gene_symbols.str.startswith("MT-")
    # Ribosomal genes
    adata.var["ribo"] = gene_symbols.str.startswith(("RPS","RPL"))
    # Hemoglobin genes
    adata.var["hb"] = gene_symbols.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt","ribo","hb"],
        inplace=True,
        log1p=False,
        percent_top=None
    )

    required_columns = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
    ]

    missing_columns = [
        column for column in required_columns
        if column not in adata.obs.columns
        ]

    if missing_columns:
        raise RuntimeError(f"QC calculation did not create expected columns: {missing_columns}")


def plot_qc_metrics(
    adata: ad.AnnData,
    fig_dir: str | Path
) -> None:
    """Violin (overall and grouped by sample) and scatter plots."""
    
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    required_metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo"
    ]

    missing_metrics = [
        metric for metric in required_metrics
        if metric not in adata.obs.columns
    ]

    if missing_metrics:
        raise ValueError(f"QC metric not found: {missing_metrics}. Run calculate_qc_metrics first.")

    # violin plots
    violins = sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False
    )
    violins.figure.savefig(
        fig_dir / "qc_violin_metrics.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(violins.figure)

    # violin plots grouped by sample
    gviolins = sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        groupby="biosample_id",
        jitter=0.4,
        rotation=90,        
        show=False,
    )
    
    if isinstance(gviolins, list):
        gviolin_figure = gviolins[0].figure
    else:
        gviolin_figure = gviolins.figure

    gviolin_figure.savefig(
        fig_dir / "qc_grouped_violin_metrics.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(gviolin_figure)

    # scatter plots
    scatter_plot = sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color=["pct_counts_mt","pct_counts_ribo"],
        show=False,
    )
    
    if isinstance(scatter_plot, list):
        scatter_figure = scatter_plot[0].figure
    else:
        scatter_figure = scatter_plot.figure

    scatter_figure.savefig(
    fig_dir / "qc_total_counts_vs_genes.png",
    dpi=300,
    bbox_inches="tight",
    )
    
    plt.close(scatter_figure)


def get_sample_qc_summary(
    adata: ad.AnnData
)-> pd.DataFrame:
    """Return sample QC metrics including: 
    
    - number of cells, 
    - median MT percentage, 
    - median number of genes, 
    - median number of counts per cell, 
    - fraction of cells over 20% MT, 
    - fraction of cells with fewer than 300 genes."""

    per_sample_qc_summary = (
        adata.obs
        .assign(
            mt_over_20=adata.obs["pct_counts_mt"] > 20,
            low_genes=adata.obs["n_genes_by_counts"] < 300
        )
        .groupby("biosample_id", observed=True)
        .agg(
            cells=("pct_counts_mt", "size"),
            median_mt=("pct_counts_mt", "median"),
            median_genes=("n_genes_by_counts", "median"),
            median_counts=("total_counts", "median"),
            frac_mt_over_20=("mt_over_20", "mean"),
            frac_low_genes=("low_genes", "mean"),
        )
        .sort_values("median_mt", ascending=False)
    )
    return per_sample_qc_summary


def mad_outlier(
    values: ArrayLike,
    nmads: float = 5.0,
    direction: Literal["lower", "upper", "both"] = "both",
    log1p: bool = False,
) -> NDArray[np.bool_]:
    """Return a Boolean mask for robust MAD-based outliers."""

    x = np.log1p(values) if log1p else np.asarray(values, dtype = float)

    median = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x-median))

    # Avoid a meaningless threshold for an invariant group
    if mad == 0 or np.isnan(mad):
        return np.zeros(len(x), dtype=bool)

    lower = median - nmads * mad
    upper = median + nmads * mad

    if direction == "lower":
        return x < lower
    if direction == "upper":
        return x > upper
    return (x < lower) | (x > upper)