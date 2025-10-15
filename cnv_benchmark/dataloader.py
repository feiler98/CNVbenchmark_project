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
from ._utility._classes import Foundation
from ._utility.dataloader_utility import _get_data_available, query_dataset, best_match
# ______________________________________________________________________________________________________________________


# Multiomics Data
# ----------------------------------------------------------------------------------------------------------------------

dict_data_section = {
                     "data_loc": str(Path(__file__).parent),
                     "requires": str(['GBC', 'RCM']),
                     "facs": "FACS"
                    }

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
    _dict_available_data = _get_data_available("data", dict_data_section)

    def __init__(self,
                 fetched_group_dict_data: dict = None,
                 selected_group: str = None):
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

========================================
         // Available Datasets
----------------------------------------""")

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
                assembly_genome = str(data_info["RCM"].stem).split(sep="__")[1]
                print(f"""        > {data_key}
           available cells | {count_cells_sample}
           assembly genome | {assembly_genome}""")

        sum_string = f"cells total     | {sum(list_count_cells)}"
        print("\n           " + "-" * len(sum_string))
        print(f"""           {sum_string}
========================================
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
        group_data = best_match(group_data, list_groups)

        # get all subsets of the group
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


# FACS Data
# ----------------------------------------------------------------------------------------------------------------------

dict_repair_facs_data = {
      "data_loc": str(Path(__file__).parent),
      "requires": str(['RCM']),
}

class FACSplus(Foundation):
    # full dictionary with all available data
    _dict_available_data = _get_data_available("facs_data", dict_repair_facs_data)

    def __init__(self,
                 df_rcm: pd.DataFrame = None,
                 df_gbc: pd.DataFrame = None,
                 facs_path_dict: dict = None):
        super().__init__()
        self.df_rcm = df_rcm
        self.df_gbc = df_gbc
        self.facs_path_dict = facs_path_dict
        self.class_obj = FACSplus
        self.error_text = """
        #################################################
        Class has not been initialized!
        Use FACSplus.fetch_matching_data(group_data: str) first.
        #################################################
        """

    # General check methods before initialization
    # -------------------------------------------
    @staticmethod
    def available_datasets(ref_genome: str = None, return_true: bool = False) -> (dict, None):
        """
        Prints a list of available datasets as overview.

        Returns
        -------
        dict
            Dictionary with all groups and their data path; can be filtered by ref_genome.
        """

        # box in the information for terminal display
        print("""

========================================
      // FACS | Available Datasets
----------------------------------------""")

        list_count_cells = []
        dict_group = {}
        for key, subdict in FACSplus._dict_available_data.items():
            # header styling
            print(f"""
    > {key}""")
            print("    " + "-" * (len(key) + 2))

            # generate the subdict information, just get the first row of the dataframe for faster loading times
            dict_accepted = {}
            for data_key, data_info in subdict.items():
                count_cells_sample = len(
                    pd.read_csv(str(data_info["RCM"]), index_col="Gene", nrows=0).columns)
                assembly_genome = str(data_info["RCM"].stem).split(sep="__")[1]
                if assembly_genome == ref_genome or ref_genome is None:
                    dict_accepted.update({data_key: data_info})
                    list_count_cells.append(count_cells_sample)
                    print(f"""        > {data_key}
               available cells | {count_cells_sample}
               assembly genome | {assembly_genome}""")
            if len(dict_accepted) == 0:
                print(f"""        > no matching data""")
            else:
                dict_group.update({key: dict_accepted})

        sum_string = f"cells total     | {sum(list_count_cells)}"
        print("\n           " + "-" * len(sum_string))
        print(f"""           {sum_string}
========================================
""")
        if return_true:
            return dict_group
        return None

    # Initialization of the facs expansion
    ####################################################################################################################
    @classmethod
    def fetch_facs_expansion(cls, rcm_data: (Path, str, pd.DataFrame), gbc_data: (Path, str, pd.DataFrame),  assembly_genome: str = None):
        """

        Parameters
        ----------
        rcm_data: PosixPath, str, pd.DataFrame
        gbc_data: PosixPath, str, pd.DataFrame
        assembly_genome: str
        """
        # check both rcm and gbc
        dict_data_check = {
            "rcm": {"data":rcm_data, "col_req":["Gene"], "df": None},
            "gbc": {"data": gbc_data, "col_req": ["CHR", "START", "END"], "df": None},
        }
        dict_accepted_facs = None


        def check_contains_list(list_query: list, list_search) -> bool:
            bool_return = True
            for items in list_query:
                if not items in list_search:
                    bool_return = False
            return bool_return

        for datatype, subdict in dict_data_check.items():
            if isinstance(subdict["data"], (Path, str)):
                import_path = Path(subdict["data"])
                if not import_path.is_file() and str(import_path).endswith((".csv", ".xlsx")):
                    raise ValueError(f"{import_path} is not a valid file! Only csv- or excel-files are accepted.")
                # check assembly_genome
                if assembly_genome is None:
                    assembly_genome = import_path.stem.split(sep="__")[1]
                df_import = pd.read_csv(import_path)
                bool_all_col_contained = check_contains_list(subdict["col_req"], list(df_import.columns))
                if not bool_all_col_contained:
                    df_import = df_import.T
                    bool_all_col_contained = check_contains_list(subdict["col_req"], list(df_import.columns))
                if not bool_all_col_contained:
                    raise ValueError(f"Required columns {subdict["col_req"]} are not in the {datatype}-file!")
                dict_data_check[datatype].update({"df":df_import.set_index(subdict["col_req"][0])})

            elif isinstance(subdict["data"], pd.DataFrame):
                # check assembly_genome
                if not isinstance(assembly_genome, str):
                    raise ValueError("""Please provide a assembly genome to the DataFrame. For orientation please look
    at the available data with FACSplus().available_datasets()""")
                df_idx_reset = subdict["data"].reset_index().drop("index", axis=1, errors="ignore")
                bool_all_col_contained = check_contains_list(subdict["col_req"], list(df_idx_reset.columns))
                if not bool_all_col_contained:
                    df_idx_reset = df_idx_reset.T
                    bool_all_col_contained = check_contains_list(subdict["col_req"], list(df_idx_reset.columns))
                if not bool_all_col_contained:
                    raise ValueError(f"Required columns {subdict["col_req"]} are not in the {datatype}-file!")
                dict_data_check[datatype].update({"df": df_idx_reset.set_index(subdict["col_req"][0])})
            else:
                raise ValueError("Datatype is not valid. Please provide a path or DataFrame.")
        dict_accepted_facs = FACSplus.available_datasets(assembly_genome, True)
        if len(dict_accepted_facs) == 0:
            raise ValueError("No match by assembly genome [file naming convention] was found for the data!")

        # returns the Multiomics DataFrame dataset and the accepted FACS-data if there were any hits by genomic assembly
        return cls(df_rcm=dict_data_check["rcm"]["df"], df_gbc=dict_data_check["gbc"]["df"], facs_path_dict=dict_accepted_facs)

    # post-initialization
    # -------------------
    def dataloader_extend_facs_hybrid(self,
                                      dataset_name: str,
                                      query_facs_data: (str, None) = None,
                                      is_facs_percent: (float, None) = None):
        """
        Expand the multiomics datasets provided by the dataloader with FACS-data.
        Integrates selected transcriptomic-data with the scverse-platform and saves the GBC- and RCM-files in a new
        folder within the dataloader directory (set in section data in the dataloader.ini).

        Parameters
        ----------
        dataset_name: str
            Name of the generated FACS multiomics hybrid dataset.
        query_facs_data: str, None
            Name of FACS data (pre-filtered by assembly genome) which fits the query by similarity score.
            If 'None', will take the whole FACS-data pool.
            Structure: group_name;data_name
                       if group_name is left empty -> all groups are searched
                       if data_name is left empty -> all datasets for group are incorporated
        is_facs_percent: float, None
            Relative percentage of the multiomics dataset as additional FACS data (randomly picked cell entities).
            Range minimum is 0 (though that would be silly to select) to max cells of the FACS data.
            1 is equal to 100% of the multiomics cell-count if available.
        """

        pass
    ####################################################################################################################

# debugging
# ---------
if __name__ == "__main__":
    pass
