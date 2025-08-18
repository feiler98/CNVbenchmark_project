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

def _get_data_available() -> dict:
    """
    Algorithm which checks and generates a dictionary of available data based on the dataloader.ini configuration-file
    section 'data'.

    Returns
    -------
    dict
        Dictionary with available multiomic datasets for every group.
    """
    path_cfg = Path(__file__).parent / "config" / "dataloader.ini"
    if not path_cfg.exists():
        print("Configuration file 'dataloader.ini' does not exist; Creating file...")
        path_cfg.touch()  # create configuration-file if it does not exist

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

    dict_available_data = _get_data_available()

    def __init__(self, fetched_group_dict_data):
        self.fetched_group_dict_data = fetched_group_dict_data

    # Initialization of the dataloader
    # --------------------------------
    @classmethod
    def fetch_data(cls, group_data: str, subset_filter: [str, list] = None):
        # list of all available group directories which meet the data criteria defined in
        # dataloader._get_data_available()
        list_groups = list(DataLoader.dict_available_data.keys())

        # check if selected group exists
        if group_data not in list_groups:
            string_groups = "\n    > ".join(list_groups)
            raise ValueError(f"""
Group is not known, please refer to the currently available groups listed below:
    > {string_groups}
            """)

        dict_subset_group = DataLoader.dict_available_data[group_data]

        # check subset_filter
        if isinstance(subset_filter, str):
            subset_filter = [subset_filter]
        elif not isinstance(subset_filter, list):
            raise ValueError("Attribute subset_filter must be the following datatypes: [str, list]")

        dict_match = {tag: dict_subset_group[tag] for tag in subset_filter if tag in list(dict_subset_group.keys())}

        if len(dict_match) == 0:
            raise ValueError(f"There are no matches for given subset_filter: {subset_filter} in {group_data}!")


    # General check methods before initialization
    # -------------------------------------------
    @staticmethod
    def available_datasets():
        # box in the information for terminal display
        print("""
        
=====================================
       // Available Datasets
-------------------------------------""")

        list_count_cells = []
        for key, subdict in DataLoader.dict_available_data.items():
            # header styling
            print(f"""
    > {key}""")
            print("    "+"-"*len(key))

            # generate the subdict information, just get the first row of the dataframe for faster loading times
            for data_key, data_info in subdict.items():
                count_cells_sample = len(pd.read_csv(str(data_info["umi_count_matrix"]), index_col="Gene", nrows=0).columns)
                list_count_cells.append(count_cells_sample)
                print(f"""        > {data_key}
           available cells | {count_cells_sample}""")

        sum_string = f"cells total     | {sum(list_count_cells)}"
        print("\n           "+"-"*len(sum_string))
        print(f"""           {sum_string}
=====================================
""")

if __name__ == "__main__":
    print(Path(__file__))
    print(_get_data_available())
    DataLoader.available_datasets()
    DataLoader.fetch_data("biaan_group")
