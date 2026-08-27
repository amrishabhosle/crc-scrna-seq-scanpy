from __future__ import annotations
from pathlib import Path
import pandas as pd
import re
import anndata as ad

GENE_ID_PATTERN = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_PATTERN = re.compile(r'gene_name "([^"]+)"')

def create_mapping(
    gtf_path: str | Path
) -> pd.DataFrame:
    """ Return unique Ensemble gene ID to symbol mapping from a GTF."""

    gtf_path = Path(gtf_path)
    if not gtf_path.exists:
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[2,8],
        names=["feature_type","attributes"],
        compression="infer",
        dtype="str"
    )
    
    genes = gtf.loc[gtf["feature_type"].eq("gene"),"attributes"].copy()
    
    # build table
    mapping = pd.DataFrame(
        {
            "ensembl_id": genes.str.extract(GENE_ID_PATTERN, expand=False),
            "gene_symbol": genes.str.extract(GENE_NAME_PATTERN, expand=False),
        }
    ).dropna(subset=["ensembl_id"])

    # remove suffix
    mapping["ensembl_id"] = mapping["ensembl_id"].str.replace(
        r"\.\d+$",
        "",
        regex=True,
    )

    mapping = mapping.drop_duplicates()
    duplicate_ids = mapping["ensembl_id"].duplicated(keep=False)
    
    if duplicate_ids.any():
        examples = mapping.loc[duplicate_ids, "ensembl_id"].head().tolist()
        raise ValueError(
            "Some Ensembl IDs map to multiple rows in the GTF. "
            f"Examples: {examples}"
        )

    return mapping.set_index("ensembl_id").sort_index()

def attach_gene_symbol(
    mapping: pd.Dataframe,
    adata: ad.AnnData
) -> ad.AnnData:
    """ Return adata with gene name annotations for Ensembl Gene IDs."""

    if not mapping.index.is_unique:
        raise ValueError("Mapping index must contain unique Ensembl IDs.")
    if not adata.var_names.is_unique:
        raise ValueError("Ensembl Gene IDs must be unique.")

    adata.var["gene_symbol"] = adata.var_names.map(mapping["gene_symbol"])

    n_matched = adata.var["gene_symbol"].notna().sum()
    n_total = adata.n_vars

    if n_matched == 0:
        raise ValueError(
            "No Ensembl IDs matched the GTF mapping. "
            "Check the genome assembly and identifier format."
        )

    adata.uns["gene_annotation"] = {
        "source_file": "Homo_sapiens.GRCh37.87.gtf.gz",
        "reference_assembly": "GRCh37",
        "ensembl_release": 87,
        "n_genes": int(n_total),
        "n_symbols_matched": int(n_matched),
        "n_symbols_unmatched": int(n_total - n_matched),
    }

    return adata