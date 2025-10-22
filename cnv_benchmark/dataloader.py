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
import scanpy as sc
from pathlib import Path
import pyomics
from pyomics import utils as ut
from ._utility._classes import Foundation
import anndata as ad
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


def _validate_data_dict(input_dict: dict, req_keys_subdict: list, add_tag: str = None) -> tuple:
    """
    Turns nested dictionary structure in 1D dictionary

    Parameters
    ----------
    input_dict: dict
        Nested dictionary of data, transformation of established data structure for dataloader.py
    req_keys_subdict: list
        List of required keys of subdicts.
    add_tag: str
        Extra string-tag for the data.

    Returns
    -------
    dict
        1D dictionary.
    """

    # validation and transformation of dataloader dict
    if add_tag is None:
        add_tag = ""
    else:
        add_tag = f"__{add_tag}"
    dict_validated = {req_key:{} for req_key in req_keys_subdict}
    set_genome = None
    for sample, subdict in input_dict.items():  # make a general function out of this and use for facs too
        if not isinstance(subdict, dict):
            raise ValueError(f"Value of dataloader_dict['{sample}'] must be opf type dict!")
        list_assembly_genome_verify = []
        for req_key in req_keys_subdict:
            if req_key not in subdict.keys():
                raise ValueError(f"{str(*req_keys_subdict)} must be in '{sample}' dict!")
            list_assembly_genome_verify.append(Path(subdict[req_key]).stem.split("__")[1])
        if len(set(list_assembly_genome_verify)) > 1:
            raise ValueError(f"Genomic Assemblies for RCM and GBC are inconsistent for {sample}!")
        if set_genome is None:
            set_genome = list_assembly_genome_verify[0]
            for req_key in req_keys_subdict:
                dict_validated[req_key].update({f"{sample}{add_tag}": subdict[req_key]})
        else:
            if set_genome == list_assembly_genome_verify[0]:
                for req_key in req_keys_subdict:
                    dict_validated[req_key].update({f"{sample}{add_tag}": subdict[req_key]})
    return dict_validated, set_genome


class FACSplus(Foundation):
    # full dictionary with all available data
    _dict_available_data = _get_data_available("facs_data", dict_repair_facs_data)

    def __init__(self,
                 mult_rcm: dict = None,
                 mult_gbc: dict = None,
                 facs_rcm: dict = None,
                 assembly_genome: str = None,
                 tag: str = None):
        super().__init__()
        self.mult_rcm = mult_rcm
        self.mult_gbc = mult_gbc
        self.facs_rcm = facs_rcm
        self.assembly_genome = assembly_genome
        self.tag = tag
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
    def fetch_facs_expansion(cls, dataloader_path_dict: dict, dataloader_add_tag: str = None):
        """

        Parameters
        ----------
        dataloader_path_dict: dict
            Dictionary of paths containing both GBC and RCM data
        dataloader_add_tag: str
            Optional tag to mark multiomics data
        """
        if dataloader_add_tag is None:
            dataloader_add_tag = "mult_data"
        dict_dataloader, set_genome = _validate_data_dict(dataloader_path_dict, ["RCM", "GBC"], dataloader_add_tag)
        dict_valid_facs = FACSplus.available_datasets(ref_genome=set_genome, return_true=True)
        dict_facs = {}
        for key, subdict in dict_valid_facs.items():
            subdict_facs, _ = _validate_data_dict(subdict, ["RCM"], key)
            if len(subdict_facs["RCM"]) > 0:
                dict_facs.update(subdict_facs["RCM"])

        return cls(mult_rcm=dict_dataloader["RCM"],
                   mult_gbc=dict_dataloader["GBC"],
                   facs_rcm=dict_facs,
                   assembly_genome=set_genome,
                   tag=dataloader_add_tag)


    # post-initialization
    # -------------------

    def get_facs_paths(self):
        # validate if initialized
        FACSplus._check_loaded(self)
        return self.facs_rcm

    def get_mult_rcm(self):
        # validate if initialized
        FACSplus._check_loaded(self)
        return self.mult_rcm

    def get_mult_gbc(self):
        # validate if initialized
        FACSplus._check_loaded(self)
        return self.mult_gbc


    def dataloader_extend_facs_hybrid(self,
                                      dataset_name: str,
                                      query_mult_data: (str, None) = None,
                                      query_facs_data: (str, None) = None,
                                      is_facs_percent: (float, int, None) = None):
        """
        Expand the multiomics datasets provided by the dataloader with FACS-data.
        Integrates selected transcriptomic-data with the scverse-platform and saves the GBC- and RCM-files in a new
        folder within the dataloader directory (set in section data in the dataloader.ini).

        Parameters
        ----------
        dataset_name: str
            Name of the generated FACS multiomics hybrid dataset.
        query_mult_data: str, None
            Name of multiomics data (pre-filtered by assembly genome) which fits the query by similarity score.
            If 'None', will take the whole multiomics-data pool.
            Space the search tags!
        query_facs_data: str, None
            Name of FACS data (pre-filtered by assembly genome) which fits the query by similarity score.
            If 'None', will take the whole FACS-data pool.
            Space the search tags!
        is_facs_percent: float, int, None
            Relative percentage of the multiomics dataset as additional FACS data (randomly picked cell entities).
            Range minimum is 0 (though that would be silly to select) to max cells of the FACS data.
            1 is equal to 100% of the multiomics cell-count if available.
        """

        # validate if initialized
        FACSplus._check_loaded(self)

        # query preparation
        list_query_mult = query_mult_data.lower().split(sep=" ")
        list_query_facs = query_facs_data.lower().split(sep=" ")

        # key list preparation
        dict_key_mult = {old.lower():old for old in list(self.mult_rcm.keys())}
        dict_key_facs = {old.lower():old for old in list(self.facs_rcm.keys())}

        # generate sets
        def get_best_match_set(query_list: list, search_list: list) -> list:
            list_out = []
            for queries in query_list:
                matches = best_match(queries, search_list, mult_match=True)
                if matches is not None:
                    list_out.extend(matches)
            return list(set(list_out))

        # multiomics-paths to adata
        # -------------------------
        if query_mult_data is not None:
            keys_mult = get_best_match_set(list_query_mult, list(dict_key_mult.keys()))
            true_mult_keys = [dict_key_mult[key] for key in keys_mult]
            paths_gbc = {k: self.mult_gbc[k] for k in true_mult_keys}
            paths_rcm = {k: self.mult_rcm[k] for k in true_mult_keys}
        else:
            paths_gbc = self.mult_gbc
            paths_rcm = self.mult_rcm


        # facs-paths to filtered adata
        # ----------------------------

        if query_facs_data is not None:
            keys_facs = get_best_match_set(list_query_facs, list(dict_key_facs.keys()))
            true_facs_keys = [dict_key_facs[key] for key in keys_facs]
            paths_facs = {k: self.facs_rcm[k] for k in true_facs_keys}
        else:
            paths_facs = self.facs_rcm

        if not isinstance(is_facs_percent, (int, float)) or is_facs_percent is not None:
            raise ValueError("Variable 'is_facs_percent' must be of the following type: int, float, None!")

        if is_facs_percent <= 0:
            is_facs_percent = 0
            adata_facs = None
        elif is_facs_percent is None:
            dict_facs_adata = {f"FACS__{tag}": sc.read_csv(path).T for tag, path in paths_facs.items()}
        else:
            list_facs_imported = []
            for _, paths in paths_facs:
                list_facs_imported.extend(list(pd.read_csv(paths, index_col="Gene", nrows=1).columns))
                list_facs_subset = ut.get_random_list_subset(list_facs_imported, len_mult_data*is_facs_percent)







    ####################################################################################################################

# debugging
# ---------
if __name__ == "__main__":
    pass
