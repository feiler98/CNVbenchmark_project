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
          "requires": ["sc_wgs_matrix", "umi_count_matrix"],
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
