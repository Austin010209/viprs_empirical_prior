from abc import ABC, abstractmethod
import numpy as np
from typing import List


class AbstractViprs(ABC):

    @abstractmethod
    def fit(self, total_epochs: int) -> List[float]:
        raise NotImplementedError("Abstract class shall not be implemented!")

    @abstractmethod
    def is_eligible(self, prior_type) -> bool:
        raise NotImplementedError("Abstract class shall not be implemented!")
