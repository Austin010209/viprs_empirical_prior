from model.abstract_viprs import AbstractViprs
from loader.chr_load_result import ChrLoadResult
from service.evaluator import Evaluator
from utils.constants import *

import numpy as np
import math
from typing import List


class EmpiricalPriorViprs(AbstractViprs):

    """ Initialize """

    def __init__(self, train_data: ChrLoadResult,
                 test_data: ChrLoadResult,
                 evaluator: Evaluator,
                 **kwargs) -> None:

        self._n_train: int = train_data.sample_size  # Sample size
        self._n_test: int = test_data.sample_size
        self._evaluator: Evaluator = evaluator

        self._m: int = train_data.n_snps  # Total SNPs
        self._k: int = train_data.n_annotations  # Total annotations
        
        # LD matrix # shape (16296, 16296)
        self._ld_train: np.ndarray = train_data.ld
        self._ld_test = test_data.ld

        self._learning_rate: float = kwargs["learning_rate"]
        assert self._learning_rate is not None

        # Marginal summary statistics of train
        self._marginal_train: np.ndarray = train_data.marginal

        # Marginal summary statistics of test
        self._marginal_test: np.ndarray = test_data.marginal  # shape (16296,)

        # M x K binary annotation matrix
        self._A: np.ndarray = train_data.annotations

        # Init parameters
        self._tau_eps: float = kwargs["tau_eps_empirical"]
        assert self._tau_eps is not None

        # var_per_snp = kwargs["var_per_snp"]
        # assert var_per_snp is not None

        # tau_init = self._k / var_per_snp
        tau_init = kwargs["tau_beta_init_empirical"]
        self._tau = np.full((self._k), tau_init, dtype=float)  # K x 1

        self._w = np.random.normal(WEIGHT_INITIALIZE_MEAN,
                                   np.sqrt(WEIGHT_INITIALIZE_VAR),
                                   size=(self._k))  # K x 1

        # Init expectations
        tau_beta_s: float = kwargs["tau_beta_s"]
        mu_beta_s: float = kwargs["mu_beta_s"]
        gamma_s: float = kwargs["gamma_s"]
        pi: float = kwargs["pi_empirical"]

        assert tau_beta_s is not None
        assert mu_beta_s is not None
        assert gamma_s is not None
        assert pi is not None

        self._tau_beta_s = np.full((self._m), tau_beta_s, dtype=float)  # M x 1
        self._mu_beta_s = np.full((self._m), mu_beta_s, dtype=float)  # M x 1
        self._gamma_s = np.full((self._m), gamma_s, dtype=float)  # M x 1
        self._pi = np.full((self._m), pi, dtype=float)  # M x 1

    """ Train CAVI model """

    def fit(self, total_epochs: int) -> List[float]:

        result = []

        for epoch in range(total_epochs):
            self._calculate_expectation()  # E Step
            self._update_parameters()  # M Step
            r_squared_test = self._evaluator.calculate_r_squared(
                self._mu_beta_s, self._marginal_test, self._ld_test)
            r_squared_train = self._evaluator.calculate_r_squared(
                self._mu_beta_s, self._marginal_train, self._ld_train)
            result.append(r_squared_test)
            print("Epoch: {}, train R-squared: {}, test R-squared: {}".format(epoch, r_squared_train, r_squared_test))

        return result

    """ E Step """

    def _calculate_expectation(self) -> None:

        for j in range(self._m):

            self._tau_beta_s[j] = self._get_tau_beta_j(j)
            self._mu_beta_s[j] = self._get_mu_beta_j(j)
            self._gamma_s[j] = self._get_gamma_j(j)
            self._pi[j] = self._get_pi_j(j)

    def _get_tau_beta_j(self, j: int) -> float:
        # Inferred precision of effect size
        r_jj = self._ld_train[j, j]
        k = np.where(self._A[j] == 1)[0]  # Changes
        tau_k_sum = 1/(np.sum(1/self._tau[k]))
        # tau_k_sum = np.sum(self._tau[k])  # Changes
        return (self._n_train * r_jj) * self._tau_eps + tau_k_sum  # Changes

    def _get_mu_beta_j(self, j: int) -> float:

        i_not_j_sum = np.sum(self._gamma_s * self._mu_beta_s *
                             self._ld_train[j, :]) - self._gamma_s[j]*self._mu_beta_s[j]*self._ld_train[j, j]
        return self._n_train*self._tau_eps/self._tau_beta_s[j] * (self._marginal_train[j] - i_not_j_sum)

    def _get_gamma_j(self, j: int) -> float:

        # Inferred PIP of effect size
        mu_j = self._get_mu_j(j)
        result = 1 / (1 + math.exp(-1 * mu_j))
        # result = 0.01 if result < 0.01 else result
        # result = 0.99 if result > 0.99 else result
        return result

    def _get_mu_j(self, j) -> float:

        tau_beta_j_s = self._tau_beta_s[j]
        mu_beta_j_s = self._mu_beta_s[j]
        pi_j = self._pi[j]  # Changes

        k = np.where(self._A[j] == 1)  # Changes
        tau_beta = np.sum(self._tau[k])  # Changes
        return math.log(pi_j / (1 - pi_j)) + 0.5 * math.log(tau_beta / tau_beta_j_s) + 0.5 * tau_beta_j_s * math.pow(mu_beta_j_s, 2)

    def _get_pi_j(self, j) -> float:

        a_j = self._A[j]  # K x 1
        x = np.sum(np.multiply(a_j, self._w))  # 1 x 1
        result = 1 / (1 + np.exp(-x))  # Changes
        # result = 0.01 if result < 0.01 else result
        # result = 0.99 if result > 0.99 else result
        return result

    """ M Step """

    def _update_parameters(self) -> None:

        # Update tau k
        for k in range(self._k):
            j = np.where(self._A[:, k] == 1)

            tau_beta_j_s = self._tau_beta_s[j]  # M_k x 1
            mu_beta_j_s = self._mu_beta_s[j]  # M_k x 1
            gamma_j_s = self._gamma_s[j]  # M_k x 1

            mu_plus_tau = np.power(mu_beta_j_s, 2) + \
                np.reciprocal(tau_beta_j_s)  # M_k x 1

            numerator = np.dot(gamma_j_s, mu_plus_tau)  # 1 x 1
            denomenator = np.sum(gamma_j_s)  # 1 x 1

            self._tau[k] = denomenator / (numerator)  # Changes

        # Update weights
        gradients = (self._gamma_s - self._pi) @ self._A  # K x 1
        self._w -= self._learning_rate * gradients  # Changes
        self._tau_eps = self.update_tau_eps()

    def update_tau_eps(self):
        u = self._mu_beta_s * self._gamma_s
        term1_last_term = (u.T @ self._ld_train @ u - np.dot(u, u)) / 2
        LD_diag = np.diag(self._ld_train)
        term1 = self._n_train / 2 \
            - float(np.dot(self._gamma_s * self._mu_beta_s, self._marginal_train*self._n_train)) \
            + 1/2*np.dot(self._gamma_s, (np.square(self._mu_beta_s) + 1/self._tau_beta_s)*LD_diag*self._n_train) \
            + term1_last_term*self._n_train
        return term1 / (self._n_train / 2)

    def is_eligible(self, prior_type) -> bool:
        return prior_type == PRIOR_TYPE_EMPIRICAL
