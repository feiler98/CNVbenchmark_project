"""
------------------------------------------------------------------------------------------------------------------------

INFERCNV_EVAL
-------------

Repository for infercnv-method results handling, correction and evaluation.

------------------------------------------------------------------------------------------------------------------------
"""

# imports
# ______________________________________________________________________________________________________________________
import pandas as pd
from pathlib import Path
# ______________________________________________________________________________________________________________________


def cnv_result_split_by_cell(path_csv: (str, Path)) -> dict:
    """
    Reduces the matrix in its redundancy (sums repeating bin values in proximity) and splits dataframe into per
    cell dataframes. Returned as dictionary.

    Parameters
    ----------
    path_csv: str | Path

    Returns
    -------
    dict
        Dictionary with cell-tag: CNV-DataFrame.
    """
    path_csv = Path(path_csv)
    if not path_csv.is_file() and path_csv.stem != ".csv":
        raise ValueError("Provided path is not a csv-file!")
    df = pd.read_csv(path_csv)
    df["CHR"] = df["CHR"].apply(lambda x: f"chr{x}")
    list_df_cells = list(df.columns)[3::]
    dict_cells_df_condensed = {}
    for cell_id in list_df_cells:
        slice_per_cell = df[["CHR", "START", "END", cell_id]]
        slice_per_cell_start = slice_per_cell.copy().drop_duplicates(subset=["CHR", cell_id], keep="first")
        slice_per_cell_end = slice_per_cell.copy().drop_duplicates(subset=["CHR", cell_id], keep="last")
        slice_per_cell_start["END"] = list(slice_per_cell_end["END"])
        slice_per_cell_start.reset_index(inplace=True, drop=True)
        dict_cells_df_condensed.update({cell_id: slice_per_cell_start})
    return dict_cells_df_condensed


if __name__ == "__main__":
    pass