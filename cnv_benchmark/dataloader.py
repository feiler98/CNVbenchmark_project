"""
------------------------------------------------------------------------------------------------------------------------

DATALOADER
----------

Handling of the Benchmarking Data for RNA InferCNV methods.
Configuration file [dataloader.ini] handles the location of the benchmarking data.
------------------------------------------------------------------------------------------------------------------------

"""

# imports
# ______________________________________________________________________________________________________________________
import pandas as pd
from pathlib import Path
import pyomics
import ast
# ______________________________________________________________________________________________________________________

def __data_available():
    path_cfg = Path(__file__).parent / "config" / "dataloader.ini"
    if not path_cfg.exists():
        print("Configuration file 'dataloader.ini' does not exist; Creating file...")
        path_cfg.touch()  # create configuration file if it does not exist

    # configuration file Object for handling the data
    cfg_obj = pyomics.GetConfig.get_config(str(path_cfg))

    # repair the data section with these standard parameters if necessary
    dict_repair_dataloader_data = {
          "data_loc": str(Path(__file__).parent),
          "requires": str(['sc_wgs_matrix', 'umi_count_matrix']),
          "overview": str(['summary', 'available_datasets']),
          "facs": "info_by_cell"
    }
    dict_data = cfg_obj.get_repair_config_section("data", dict_repair_dataloader_data)

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
    # priority has if the DNA and RNA files are present, all other information regarding overview and facs data are optional
    dict_data_overview = {}
    requires_list =  ast.literal_eval(dict_data["requires"])
    for data_name, p in dict_data_dir.items():
        dict_available_genomic_files = {p_csv.stem: p_csv for p_csv in p.glob("*.csv")} # get list of all csv file paths --> standard format of count matrices
        set_data_tags = set([p.stem.split(sep="__")[-1] for p in list(dict_available_genomic_files.values())])  # get unique dataset identifier tag
        dict_accepted = {}
        for tags in set_data_tags:
            flag_keep = True  # keep the dataset or remove it
            dict_paths_files = {}
            for data_type in requires_list:  # must contain all required data files
                if f"{data_type}__{tags}" in list(dict_available_genomic_files.keys()):
                    dict_paths_files[data_type] = dict_available_genomic_files[f"{data_type}__{tags}"]
                else:
                    flag_keep = False
            if flag_keep:
                dict_accepted[tags] = dict_paths_files
        if len(dict_accepted) > 0:  # in order to keep a "lab group", must at least contain 1 dataset
            dict_data_overview[data_name] = dict_accepted
    return dict_data_overview





class DataLoader:

    def __init__(self):
        pass

    # getting the necessary config information


    @classmethod
    def fetch_data(cls):
        pass

    @staticmethod
    def available_datasets():
        pass


if __name__ == "__main__":
    print(Path(__file__))
    print(__data_available())

