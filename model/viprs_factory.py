from model.abstract_viprs import AbstractViprs
from model.empirical_prior_viprs import EmpiricalPriorViprs
from model.sas_prior_viprs import SasPriorViprs
from loader.chr_load_result import ChrLoadResult
from service.evaluator import Evaluator


class ViprsFactory:

    """ Initialize """

    def __init__(self, train_data: ChrLoadResult,
                 test_data: ChrLoadResult,
                 evaluator: Evaluator,
                 **kwargs) -> None:
        self._viprs = []
        children = AbstractViprs.__subclasses__()
        if len(children) > 0:
            for child in children:
                self._viprs.append(
                    child(train_data, test_data, evaluator, **kwargs))

    """ Public methods """

    def get_viprs(self, prior_type: str) -> AbstractViprs:

        if len(self._viprs) == 0:
            return None
        for viprs in self._viprs:
            if viprs.is_eligible(prior_type):
                return viprs

        return None
