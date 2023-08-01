from loader.chr_data_loader import ChrLoadResult

import numpy as np
import magenpy as mgp
from utils.constants import *


class ChrDataProcessor:

    """ Intialize """

    def __init__(self, chr: int, fold: int, is_train=True) -> None:

        self._chr = chr

        sumstats_files = TRAIN_SUMSTATS_FILES if is_train else TEST_SUMSTATS_FILES
        sumstats_format = SUMSTATS_FORMAT_PLINK if is_train else SUMSTATS_FORMAT_MAGENPY

        self._data_loader = mgp.GWADataLoader(ld_store_files=LD_STORE_FILES.format(chr),
                                              sumstats_files=sumstats_files.format(
                                                  fold),
                                              sumstats_format=sumstats_format,
                                              annotation_files=ANNOTATION_FILES.format(
                                                  chr),
                                              annotation_format=ANNOTATION_FORMAT_LDSC)

    """ Public methods """

    def get_snps(self) -> np.ndarray:
        return self._data_loader.snps[self._chr]

    def filter(self, common_ids: np.ndarray) -> ChrLoadResult:

        # Filter by common SNPs
        sample_size = self._data_loader.sample_size
        n_snps = common_ids.shape[0]
        n_annotations = self._data_loader.annotation[self._chr].n_annotations

        ld_csr = self._data_loader.get_ld_matrices()[self._chr].to_csr_matrix()
        ld = ld_csr.toarray()[common_ids][:, common_ids]

        marginal = self._data_loader.sumstats_table[self._chr].marginal_beta
        marginal = marginal[common_ids]

        # standardizing marginal
        marginal_mean = np.mean(marginal, axis = 0)
        marginal_std = np.std(marginal, axis = 0)
        marginal = (marginal - marginal_mean) / marginal_std
        marginal *= np.sqrt(1/sample_size)

        annotations = self._data_loader.annotation[self._chr].values()
        annotations = annotations[common_ids]

        # Keep binary annotations
        annot_masks = np.where((annotations != 0.0) & (annotations != 1.0))[1]
        keep_annot_cols = [i for i in range(
            annotations.shape[1]) if i not in annot_masks]
        keep_annotations = annotations[:, keep_annot_cols]
        keep_n_annotations = len(keep_annot_cols)

        assert np.all(np.in1d(keep_annotations, [0., 1.]))

        print("Sample size: {}, total SNPs: {}, total annotations: {}, LD: {}, marginal: {}, annotations: {}".format(
            sample_size, n_snps, keep_n_annotations, ld.shape, marginal.shape, keep_annotations.shape))

        return ChrLoadResult(sample_size=sample_size,
                             n_snps=n_snps,
                             n_annotations=keep_n_annotations,
                             ld=ld,
                             marginal=marginal,
                             annotations=keep_annotations)
