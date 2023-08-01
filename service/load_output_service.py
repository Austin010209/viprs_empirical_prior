from utils.constants import *
from typing import Dict, List

import numpy as np
import os
import re


class LoadOutputService:

    def load(self, output_files: List[str]) -> Dict[str, Dict[str, List[np.ndarray]]]:

        result = None

        if output_files is not None and len(output_files) > 0:

            result: Dict[str, Dict[str, List[np.ndarray]]] = dict()

            result["21"] = {"True": [],
                            "False": []}
            result["22"] = {"True": [],
                            "False": []}

            for file_name in output_files:

                if not file_name.endswith(".csv"):
                    continue

                chr = re.search(
                    r'chr_(.*?)_use_annot', file_name).group(1)
                use_annot = re.search(
                    r'use_annot_(.*?).csv', file_name).group(1)

                # Load file
                file_path = os.path.join(OUTPUT_DIR, file_name)
                r_squared_result = np.genfromtxt(file_path)

                result_list = result[chr][use_annot]
                result_list.append(r_squared_result)

        return result
