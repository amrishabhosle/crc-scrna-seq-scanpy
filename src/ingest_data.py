from __future__ import annotations
from pathlib import Path
from scipy.io import mmread
import pandas as pd
import numpy as np
import anndata as ad

def load_expression_atlas_raw_counts(
        data_dir: str | Path,
        accession: str
) -> ad.AnnData:
    """
    Loads an EBI Expression Atlas sc raw-count Matrix Market bundle into AnnData.

    Parameters
    ----------
    data_dir
        Folder that contains the downloaded EBI files.
    accession
        Study accession used in the EBI file names, such as "E-MTAB-8410".
    metadata_filename
        Optional metadata file name. If omitted, use the standard EBI name.

    Returns
    -------
    ad.AnnData
        An AnnData object with cells as rows and genes as columns.
    """
    data_dir = Path(data_dir)

    # Build paths
    matrix_path = data_dir / f"{accession}.aggregated_counts.mtx.gz"
    rows_path = data_dir / f"{accession}.aggregated_counts.mtx_rows.gz"
    cols_path = data_dir / f"{accession}.aggregated_counts.mtx_cols.gz"
    metadata_path = data_dir / f"{accession}.cell_metadata.tsv"

    # Store paths that must exist before loading begins.
    required_paths = [matrix_path, rows_path, cols_path, metadata_path]

    # Stop if files not found
    missing_paths = [str(path) for path in required_paths if not path.exists()]
    if missing_paths:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing_paths))

    # Read the Matrix Market file
    # Source orientation is genes x cells
    # tocsr() converts it to compressed sparse row format for efficient use.
    genes_by_cells = mmread(matrix_path).tocsr()

    # Read gene ID file
    gene_table = pd.read_csv(
        rows_path,
        sep=r"\s+",
        header=None,
        names=["ensembl_id", "source_gene_id"],
        dtype=str,
    )

    # Read cell ID file
    cell_ids = pd.read_csv(
        cols_path,
        sep="\t",
        header=None,
        names=["cell_id"],
        dtype=str,
    )["cell_id"]

    # Read metadata
    metadata = pd.read_csv(
        metadata_path,
        sep="\t",
        index_col="id",
        dtype=str,
    )

    # Confirm the matrix dimensions match the gene and cell ID lists.
    # A mismatch means the files should not be combined.
    if genes_by_cells.shape != (len(gene_table), len(cell_ids)):
        raise ValueError(
            "Matrix dimensions do not match the supplied gene and cell ID files. "
            f"Matrix: {genes_by_cells.shape}; "
            f"genes: {len(gene_table)}; cells: {len(cell_ids)}."
        )

    # Require a unique Ensembl ID for every gene column in AnnData.
    if not gene_table["ensembl_id"].is_unique:
        raise ValueError("Ensembl gene IDs are not unique.")

    # Require a unique barcode or cell identifier for each matrix column.
    if not cell_ids.is_unique:
        raise ValueError("Cell IDs are not unique.")

    # Require one unique metadata record per cell ID.
    if not metadata.index.is_unique:
        raise ValueError("Metadata cell IDs are not unique.")

    # Identify matrix cell IDs that cannot be found in metadata.
    missing_metadata = cell_ids[~cell_ids.isin(metadata.index)]

    if not missing_metadata.empty:
        raise ValueError(
            f"{len(missing_metadata)} matrix cell IDs have no matching metadata."
        )

    # Construct the AnnData object
    adata = ad.AnnData(
        X = genes_by_cells.T.tocsr(),
        obs = metadata.loc[cell_ids].copy(),
        var = gene_table.set_index("ensembl_id").copy(),
    )

    # Explicitly set the observation names to the ordered matrix cell IDs.
    adata.obs_names = cell_ids.to_numpy()

    # Explicitly set the variable names to the ordered Ensembl IDs.
    adata.var_names = gene_table["ensembl_id"].to_numpy()

    # Check values
    counts = adata.X.data

    # Raw counts cannot be negative.
    if np.any(counts < 0):
        raise ValueError("Expression matrix contains negative values.")

    fraction_of_fractions = float(np.mean(~np.isclose(counts, np.round(counts))))

    # Store human-readable provenance inside the resulting H5AD file.
    adata.uns["data_provenance"] = {
        "accession": accession,
        "source": "EMBL-EBI Single Cell Expression Atlas",
        "matrix_file": matrix_path.name,
        "matrix_type": "unfiltered raw counts",
        "source_orientation": "genes_by_cells",
        "anndata_orientation": "cells_by_genes",
        "fraction_non_integer_nonzero_values": fraction_of_fractions,
        "gene_identifier": "Ensembl gene ID",
    }

    # Create a counts layer for raw counts
    adata.layers["counts"] = adata.X.copy()

    return adata

def write_h5ad(
        adata: ad.AnnData,
        output_path: str | Path
) -> None:

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")