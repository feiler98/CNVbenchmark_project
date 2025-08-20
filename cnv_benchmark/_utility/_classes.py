"""
------------------------------------------------------------------------------------------------------------------------

_CLASSES
--------

Base Classes for INHERITANCE
------------------------------------------------------------------------------------------------------------------------

"""

# imports
# ______________________________________________________________________________________________________________________
import inspect
# ______________________________________________________________________________________________________________________

class Foundation:

    def __init__(self):
        self.class_obj = Foundation
        self.error_text = """
#################################################
    Some example text here for the ValueError
#################################################
"""

    # help to get a full overview of what the class has to offer
    def help(self) -> None:
        # rename words for interpretability’s sake
        dict_replace = {
            "bound": "classmethod",
            "function": "method"
        }
        method_list = inspect.getmembers(self.class_obj)
        print("""
------------------------------
| Available callable methods |
------------------------------""")
        for items in method_list:
            if not items[0].startswith('_'):
                string_description = str(items[1]).replace("<", "").split(" ")[0]
                if string_description in dict_replace.keys():
                    string_description = dict_replace[string_description]
                print(
                    f"    > {str(self.class_obj).replace("'>", "").split(".")[-1]}().{items[0]}() -> {string_description}")


    # check if the Class has been initialized
    # ---------------------------------------
    def _check_loaded(self):
        for keys, values in self.__dict__.items():
            if values is None:
                raise ValueError(self.error_text)