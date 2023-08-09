import re
from collections import defaultdict
import pandas as pd


def parse_sdf(file_name: str):
    # Only tested on Chebi, ignores atom position matrix

    with open(file_name, "r") as file:
        lines = file.readlines()

    sdf_id_pattern = re.compile("^> <(.*?)>")

    value_dict_list = list()
    values_dict = defaultdict(list)
    var_name = ""

    for index in range(len(lines) - 1):
        if lines[index].startswith("$$$$"):
            value_dict_list.append(values_dict)
            values_dict = defaultdict(list)
            var_name = ""
        elif lines[index].startswith("> "):
            var_name = re.search(sdf_id_pattern, lines[index]).group(1)
        elif lines[index].strip() and var_name:
            value = lines[index].strip()
            values_dict[var_name].append(value)
    value_dict_list.append(values_dict)

    df_chebi_sdf = pd.DataFrame.from_dict(value_dict_list).applymap(
        lambda el: el[0] if (isinstance(el, list) and len(el) == 1) else el
    )
    return df_chebi_sdf
