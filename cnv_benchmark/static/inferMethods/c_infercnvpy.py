"""
------------------------------------------------------------------------------------------------------------------------

C_INFERCNVPY
------------

Code for scverse infercnvpy
Documentation: https://infercnvpy.readthedocs.io/en/latest/

------------------------------------------------------------------------------------------------------------------------
"""

# imports
# ______________________________________________________________________________________________________________________
import pandas as pd
from pathlib import Path
from pyomics.utils import benchmark_method
import infercnvpy as cnv
from anndata import AnnData
# ______________________________________________________________________________________________________________________


def construct_adata_from_csv(path_csv: Path, kwargs: dict) -> AnnData:
    """
    Builds an AnnData-object from an imported .csv RNA UMI-matrix.

    Parameters:
    -----------
    path_csv: pathlib.Path
        Posixpath provided by DataLoader.
    kwargs:
        Key-word-arguments (= kwargs) for infercnvpy.tl.infercnv() function.

    Returns
    -------
    AnnData
        Transformed .csv file for downstream CSV prediction by infercnvpy.
    """
    pass


def run_py_infercnv(adata: AnnData, kwargs: dict = dict) -> pd.DataFrame:
    """
    Main function for running infercnvpy for benchmarking.

    Parameters
    ----------
    adata: AnnData
        AnnData object with the RNA UMI-matrix.
    kwargs: dict
        Key-word-arguments (= kwargs) for infercnvpy.tl.infercnv() function.

    Returns
    -------
    pd.DataFrame
        Returns inferred copy number variations as table.
    """

    cnv.tl.infercnv(adata, **kwargs)
    cnv_X = adata.obsm["X_cnv"].toarray()

    # shape of the numpy-array cnv_X
    shape_tuple = cnv_X.shape
    cnv_idx = list(adata.obs.index)

    # col_tags
    chr_start_pos_list = list(adata.uns["cnv"]["chr_pos"].values())
    chr_tags = list(adata.uns["cnv"]["chr_pos"].keys())
    chr_start_pos_list_pop_first = chr_start_pos_list[1::]
    chr_start_pos_list_pop_first.append(shape_tuple[1])
    dict_chr_end_start_pos = {chr_tag :range(int(e ) -int(s)) for chr_tag, s, e in zip(chr_tags, chr_start_pos_list, chr_start_pos_list_pop_first)}

    cnv_col = []
    for chr_tag, range_list in dict_chr_end_start_pos.items():
        for pos in range_list:
            cnv_col.append(f"{chr_tag}__rel_pos_{pos}")
    return pd.DataFrame(data=cnv_X, columns=cnv_col, index=cnv_idx).T