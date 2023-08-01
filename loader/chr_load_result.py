import numpy as np


class ChrLoadResult:

    """ Initialize """

    def __init__(self, sample_size: int,
                 n_snps: int,
                 n_annotations: int,
                 ld: np.ndarray,
                 marginal: np.ndarray,
                 annotations: np.ndarray) -> None:

        self._sample_size = sample_size
        self._n_snps = n_snps
        self._n_annotations = n_annotations
        self._ld = ld
        self._marginal = marginal
        self._annotations = annotations

    """ Getters """

    @property
    def sample_size(self) -> int:
        return self._sample_size

    @property
    def n_snps(self) -> int:
        return self._n_snps

    @property
    def n_annotations(self) -> int:
        return self._n_annotations

    @property
    def ld(self) -> np.ndarray:
        return self._ld

    @property
    def marginal(self) -> np.ndarray:
        return self._marginal

    @property
    def annotations(self) -> np.ndarray:
        return self._annotations
