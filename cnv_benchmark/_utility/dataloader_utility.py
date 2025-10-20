"""
------------------------------------------------------------------------------------------------------------------------

DATALOADER_UTILITY
------------------

Utility for the multiomics and facs data.
General functions which optimize handling of the data

------------------------------------------------------------------------------------------------------------------------
"""

# imports
# ______________________________________________________________________________________________________________________
import ast
import re
from pathlib import Path
from Bio import Align
import pyomics
# ______________________________________________________________________________________________________________________


def _get_data_available(section: str, dict_repair_dataloader_data: dict) -> dict:
    """
    Algorithm which checks and generates a dictionary of available data based on the dataloader.ini configuration-file
    section 'data'.

    Parameters
    ----------
    section: str
        Name of the target section within the config-file.
    dict_repair_dataloader_data: dict
        Standard parameter dictionary of the config-file section.

    Returns
    -------
    dict
        Dictionary with available multiomic datasets for every group.
    """

    path_cfg = Path(__file__).parent.parent / "config" / "dataloader.ini"
    if not path_cfg.exists():
        print("Configuration file 'dataloader.ini' does not exist; Creating file...")
        path_cfg.touch()  # create configuration-file if it does not exist

    # configuration file Object for handling the data
    cfg_obj = pyomics.GetConfig.get_config(str(path_cfg))

    # repair the data section with these standard parameters if necessary
    dict_data = cfg_obj.get_repair_config_section(section, dict_repair_dataloader_data)

    # check and collect the available information about the config-set directory
    path_data = Path(dict_data["data_loc"]).resolve()
    # Error 1: directory does not exist
    if not path_data.exists():
        raise ValueError(f"Path {path_data} does not exist!")

    # get a dictionary with the directory names and paths --> the directory name will be the listed dataset name!
    dict_data_dir = {p.name: p for p in path_data.iterdir() if p.is_dir()}

    # Error 2: no dir in parent dir
    if len(dict_data_dir) < 1:
        raise ValueError(f"{path_data} does not contain any directories. Please recheck that the set 'data_loc' path contains")

    # structure for checking the available datasets:
    # ----------------------------------------------
    # priority: the DNA and RNA files are present, all other information regarding overview and facs data are optional
    dict_data_overview = {}
    requires_list =  ast.literal_eval(dict_data["requires"])
    for data_name, p in dict_data_dir.items():
        dict_available_genomic_files = {p_csv.stem: p_csv for p_csv in p.glob("*.csv")} # get list of all csv file paths --> standard format of count matrices
        set_data_tags = set([p.stem.split(sep="__")[0] for p in list(dict_available_genomic_files.values())])  # get unique dataset identifier tag
        dict_accepted = {}
        for tags in set_data_tags:
            flag_keep = True  # keep the dataset or remove it
            dict_paths_files = {}
            for data_type in requires_list:  # must contain all required data files

                # Regular Expression for taking the chromosome into account; .group() for getting result as string!
                list_keys = list(dict_available_genomic_files.keys())
                re_tag = f"{tags}__[a-zA-Z]+_[0-9]+__{data_type}"
                list_match = [re.match(re_tag, item).group() for item in list_keys if re.match(re_tag, item) is not None]
                if len(list_match) > 0:
                    dict_paths_files[data_type] = dict_available_genomic_files[list_match[0]]
                else:
                    flag_keep = False
            if flag_keep:
                dict_accepted[tags] = dict_paths_files
        if len(dict_accepted) > 0:  # in order to keep a "lab group", must at least contain 1 dataset
            dict_data_overview[data_name] = dict_accepted
    return dict_data_overview


def best_match(query: str, search_list: list, mult_match: bool = False) -> (str, list, None):
    """
    Matching a query to a list of strings. Returns the best match.
    Raises an error if too ambigious and multiple have the same score.

    Parameters
    ----------
    query: str
        The query the list of sequences should be checked for maximized match.
    search_list: list
        List of strings which are searched.
    mult_match: bool
        Enables multiple matches; returns all items with max score.
    Returns
    -------
    str, list, None
        The best match or multiple matches if mult_match is True
    """

    aligner = Align.PairwiseAligner(match_score=1.0)
    aligner.mode = "local"
    aligner.open_gap_score = -1
    aligner.extend_gap_score = 0
    list_scores = [aligner.score(str(target), query) for target in search_list]
    if not mult_match:
        if list_scores.count(max(list_scores)) > 1:
            raise ValueError("Match of query is too ambigious. Please be more precise with your query!")
        hit = search_list[list_scores.index(max(list_scores))]
    else:
        if max(list_scores) > 0:
            scores_idx = [i for i, x in enumerate(list_scores) if x == max(list_scores)]
            hit = [search_list[idx] for idx in scores_idx]
        else:
            hit = None
    return hit


def query_dataset(query: str, data_path_dict: dict) -> Path:
    """
    Finds the best file-path by query. Case insensitive!

    Parameters
    ----------
    query: str
        The query the list of sequences should be checked for maximized match.
    data_path_dict: dict
        Nested dictionary structure, special implementation for the dataloader.py script.

    Returns
    -------
    PosixPath
        Technically it returns the value of the matching nested dictionary.
    """

    lquery = query.lower()
    dict_lqueried = {}
    for key, subdict in data_path_dict.items():
        for subkey, path in subdict.items():
            dict_lqueried.update({f"{key} {subkey}".lower():path})
    best_match_key = best_match(lquery, list(dict_lqueried.keys()))
    # return best matching file-path to the query
    return dict_lqueried[best_match_key]["RCM"]


# debugging
if __name__ == "__main__":
    pass










