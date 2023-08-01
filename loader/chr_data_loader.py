from loader.chr_load_result import ChrLoadResult
from loader.chr_data_processor import ChrDataProcessor
from utils.constants import *

from typing import Dict, Tuple
import numpy as np


class ChrDataLoader:

    """ Initialize """

    # def __init__(self) -> None:

    #     # Key: (chromosome, fold)
    #     # Value: (train_data, test_data)
    #     self._result_dict: Dict[Tuple[int, int],
    #                             Tuple[ChrLoadResult, ChrLoadResult]] = dict()

    """ Public method """

    def load(self, chr: int, fold: int) -> Tuple[ChrLoadResult, ChrLoadResult]:

        print("Load chr: {}, fold: {}".format(chr, fold))

        # Early return if load result exists
        # result = self._result_dict.get((chr, fold))
        # if result is not None:
        #     return result

        # Load data
        train_processor = ChrDataProcessor(chr, fold, True)
        train_snps = train_processor.get_snps()

        test_processor = ChrDataProcessor(chr, fold, False)
        test_snps = test_processor.get_snps()

        print("Train: {}, shape: {}".format(train_snps, train_snps.shape))
        print("Test: {}, shape: {}".format(test_snps, test_snps.shape))

        # Find common SNPs
        common_values = np.unique(np.intersect1d(train_snps, test_snps))
        print("Found common SNPs: {}, shape: {}".format(
            common_values, common_values.shape))

        train_snps_ids = np.where(np.in1d(train_snps, common_values))[0]
        test_snps_ids = np.where(np.in1d(test_snps, common_values))[0]

        # Filter by common SNPs
        train_result = train_processor.filter(train_snps_ids)
        test_result = test_processor.filter(test_snps_ids)

        # self._result_dict[(chr, fold)] = (train_result, test_result)

        return train_result, test_result
