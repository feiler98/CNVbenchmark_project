from pyomics.utils import get_project_dir
from .dataloader import *

# callable information about cnv_benchmark
########################################################################################################################
import importlib.metadata

def __get_txt_path(txt_name):
    path_main = get_project_dir(Path(__file__).parent.resolve(), "CNVbenchmark_project")
    txt_path = path_main / "cnv_benchmark" / "static" / txt_name
    return txt_path

def __print_txt_line_by_line(path: str):
    """
    Proper Display of information callables in e.g. Jupyter Notebooks

    Paramters
    ---------
    path: str
        Absolute string path of .txt file.

    Returns
    -------
        Prints line by line of .txt file into the terminal.
    """
    with open(path, "r") as f:  # r for reading
        for lines in f:
            print(lines.replace("\n", ""))  # .strip() removes empty spaces at the beginning of the line

__version__ = f"cnv_benchmark {importlib.metadata.version("cnv_benchmark")}"


def quickstart():
    __print_txt_line_by_line(__get_txt_path("quickstart.txt"))

