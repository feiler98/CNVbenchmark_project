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
from ._utility._classes import Foundation
import re
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
          "requires": str(['GBC', 'RCM']),
          "facs": "FACS"
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


class DataLoader(Foundation):
    """
    Handling of the data access for the benchmarking. Dependent on the given dataloader.ini configuration file.

    Attributes
    ----------
    fetched_group_dict_data: dict
        Dictionary of the available datasets, including the paths for the required datafiles.
    selected_group: str
        Tag-name of group the selected data originates from.
    """

    # full dictionary with all available data
    _dict_available_data = _get_data_available()

    def __init__(self, fetched_group_dict_data: dict = None, selected_group: str = None):
        super().__init__()
        self.fetched_group_dict_data = fetched_group_dict_data
        self.selected_group = selected_group
        self.class_obj = DataLoader
        self.error_text = """
        #################################################
        Class has not been initialized!
        Use Dataloader.fetch_data(group_data: str) first.
        #################################################
        """

    # General check methods before initialization
    # -------------------------------------------
    @staticmethod
    def available_datasets() -> None:
        """
        Prints a list of available datasets as overview.

        Retruns
        -------
        None
        """

        # box in the information for terminal display
        print("""

=====================================
       // Available Datasets
-------------------------------------""")

        list_count_cells = []
        for key, subdict in DataLoader._dict_available_data.items():
            # header styling
            print(f"""
    > {key}""")
            print("    " + "-" * (len(key)+2))

            # generate the subdict information, just get the first row of the dataframe for faster loading times
            for data_key, data_info in subdict.items():
                count_cells_sample = len(
                    pd.read_csv(str(data_info["RCM"]), index_col="Gene", nrows=0).columns)
                list_count_cells.append(count_cells_sample)
                print(f"""        > {data_key}
           available cells | {count_cells_sample}""")

        sum_string = f"cells total     | {sum(list_count_cells)}"
        print("\n           " + "-" * len(sum_string))
        print(f"""           {sum_string}
=====================================
""")


    # Initialization of the dataloader
    ####################################################################################################################
    @classmethod
    def fetch_data(cls, group_data: str, subset_filter: (str, list) = None):
        """
        The @classmethod; initializes the class by specifying a specific group and their data.
        It is recommended to first check the DataLoader().available_datasets() before initializing the class.

        Parameters
        ----------
        group_data: str
            Name tag of the group.
        subset_filter: (str, list)
            Either a name of a specific dataset or a list of dataset names.

        Returns
        -------
        None
            Post-initialization methods are used to access the data:
            get_as_dataframe(); get_as_path(), get_group_info(), get_group_summary()
        """

        # list of all available group directories which meet the data criteria defined in
        # dataloader._get_data_available()
        list_groups = list(DataLoader._dict_available_data.keys())

        # check if selected group exists
        if group_data not in list_groups:
            string_groups = "\n    > ".join(list_groups)
            raise ValueError(f"""
Group is not known, please refer to the currently available groups listed below:
    > {string_groups}
            """)

        dict_subset_group = DataLoader._dict_available_data[group_data]

        # check subset_filter
        if isinstance(subset_filter, str):
            subset_filter = [subset_filter]
        elif not isinstance(subset_filter, (str, list)) and subset_filter is not None:
            raise ValueError("Attribute subset_filter must be the following datatypes: (str, list)")

        if subset_filter is None:
            dict_match = dict_subset_group
        else:
            dict_match = {tag: dict_subset_group[tag] for tag in subset_filter if tag in list(dict_subset_group.keys())}

        if len(dict_match) == 0:
            raise ValueError(f"There are no matches for given subset_filter: {subset_filter} in {group_data}!")

        return cls(fetched_group_dict_data=dict_match, selected_group=group_data)

    # post-initialization
    # -------------------
    def get_as_dataframe(self) -> dict:
        """
        Get the data as pd.DataFrame.

        Returns
        -------
        dict
            Dictionary of sorted pd.DataFrames.
        """

        Foundation._check_loaded(self)
        for key, data in self.fetched_group_dict_data.items():
            self.fetched_group_dict_data[key]['GBC'] = pd.read_csv(str(self.fetched_group_dict_data[key]['GBC']), index_col="CHR")
            self.fetched_group_dict_data[key]['RCM'] = pd.read_csv(str(self.fetched_group_dict_data[key]['RCM']), index_col="Gene")
        return self.fetched_group_dict_data

    def get_as_path(self) -> dict:
        """
        Get the pathlib.Path object for each data-file.

        Returns
        -------
        dict
            Dictionary of sorted file-paths for low memory usage and individualized import.
        """

        Foundation._check_loaded(self)
        return self.fetched_group_dict_data

    def _group_dir_path(self) -> Path:
        path_cfg = Path(__file__).parent / "config" / "dataloader.ini"
        # configuration file Object for handling the data
        cfg_obj = pyomics.GetConfig.get_config(str(path_cfg))
        path_group = Path(cfg_obj.return_section("data")["data_loc"]) / self.selected_group
        return path_group

    def get_group_info(self) -> (pd.DataFrame, None):
        """
        Information about the group/authors of the datasets.

        Returns
        -------
        pd.DataFrame
            Table with information of the authors and the link to the data-associated publication.
        """
        Foundation._check_loaded(self)
        path_group = DataLoader._group_dir_path(self) / f"{self.selected_group.replace("_group", "")}_summary.xlsx"
        if path_group.exists():
            return pd.read_excel(str(path_group), index_col="index")
        else:
            print(f"Summary file for {self.selected_group} is not available.")
            return None

    def get_data_summary(self) -> (pd.DataFrame, None):
        """
        Information about the loaded datasets; is dependen on the chosen subset_filter of the @classmethod.

        Returns
        -------
        pd.DataFrame
            Overview table with information about the individual datasets.
        """

        Foundation._check_loaded(self)
        path_group = DataLoader._group_dir_path(self) / f"{self.selected_group.replace("_group", "")}_availableDatasets.xlsx"
        if path_group.exists():
            selected_datasets = list(self.fetched_group_dict_data.keys())
            df = pd.read_excel(str(path_group), index_col="index")
            filtered_cols = [c for c in selected_datasets if c in list(df.columns)]
            return df[filtered_cols] if len(filtered_cols) > 0 else None
        else:
            print(f"Available_dataset file for {self.selected_group} is not available.")
            return None
    ####################################################################################################################


# debugging
# ---------
if __name__ == "__main__":
    pass
    # print(Path(__file__))
    # print(_get_data_available())
    # DataLoader.available_datasets()
    # DataLoader().help()
    # wu_dataloader = DataLoader.fetch_data("wu_group", subset_filter="GBM")
    # print(wu_dataloader.get_data_summary())
    # print(wu_dataloader.get_group_info())
    # DataLoader().available_datasets()

    # must raise ValueError!
    # DataLoader().get_data_summary()
