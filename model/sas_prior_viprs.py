from model.abstract_viprs import AbstractViprs
from loader.chr_load_result import ChrLoadResult
from service.evaluator import Evaluator
from utils.constants import *

from scipy.special import expit
import numpy as np


class SasPriorViprs(AbstractViprs):

    def __init__(self, train_data: ChrLoadResult,
                 test_data: ChrLoadResult,
                 evaluator: Evaluator,
                 **kwargs):

        self.M = train_data.n_snps
        self.N_train = train_data.sample_size
        self.N_test = test_data.sample_size
        self.LD = train_data.ld
        self._test_ld = test_data.ld
        self.marginal_stats = train_data.marginal
        self._marginal_test = test_data.marginal
        self._evaluator = evaluator

        self.mu_star = np.full((self.M), kwargs["mu_beta_s"], dtype=float)
        self.tao_beta_star = np.full(
            (self.M), kwargs["tau_beta_s"], dtype=float)
        self.gamma_star = np.full((self.M), kwargs["gamma_s"], dtype=float)
        self.tao_epsilon = kwargs["tau_eps_sas"]
        self.tao_beta = kwargs["tau_beta"]
        self.pi = kwargs["pi_sas"]

    def E_step(self):
        for j in range(self.M):
            self.tao_beta_star[j] = self.N_train * \
                self.LD[j][j] * self.tao_epsilon + self.tao_beta
            i_not_j_sum = np.sum(self.gamma_star * self.mu_star *
                                 self.LD[j, :]) - self.gamma_star[j]*self.mu_star[j]*self.LD[j, j]
            self.mu_star[j] = self.N_train*self.tao_epsilon / \
                self.tao_beta_star[j] * (self.marginal_stats[j] - i_not_j_sum)
            u_j = np.log(self.pi/(1-self.pi)) + 0.5*np.log(self.tao_beta /
                                                           self.tao_beta_star[j]) + 0.5*self.tao_beta_star[j]*(self.mu_star[j])**2
            self.gamma_star[j] = expit(u_j)
            self.gamma_star[j] = max(min(0.99, self.gamma_star[j]), 0.01)

    def M_step(self):
        tao_beta_inv = np.dot(self.gamma_star, np.square(
            self.mu_star) + 1/self.tao_beta_star) / np.sum(self.gamma_star)
        self.tao_beta = 1/tao_beta_inv
        self.pi = np.sum(self.gamma_star) / self.M

    def evaluate_ELBO(self):

        term1_last_term = 0
        for j in range(self.M):
            for k in range(j+1, self.M):
                term1_last_term += self.gamma_star[j]*self.mu_star[j] * \
                    self.gamma_star[k]*self.mu_star[k] * \
                    self.LD[k][j]*self.N_train

        LD_diag = np.diag(self.LD)
        term1 = self.N_train/2*np.log(self.tao_epsilon) - self.tao_epsilon/2*self.N_train + \
            float(self.tao_epsilon*np.dot(self.gamma_star * self.mu_star, self.marginal_stats*self.N_train)) - \
            self.tao_epsilon/2*np.dot(self.gamma_star, (np.square(self.mu_star) + 1 /
                                      self.tao_beta_star)*LD_diag*self.N_train) - self.tao_epsilon*term1_last_term

        term2 = -0.5*self.M*np.log(2*np.pi/self.tao_beta) - self.tao_beta/2 * np.dot(
            self.gamma_star, np.square(self.mu_star) + 1/self.tao_beta_star)
        term3 = -0.5*self.M*np.log(2*np.pi/self.tao_beta) - \
            0.5*np.log(self.tao_beta)*np.sum(self.gamma_star)
        term4 = np.log(self.pi)*np.sum(self.gamma_star) + \
            np.log(1-self.pi)*np.sum(1-self.gamma_star)
        term5 = np.sum(self.gamma_star * np.log(self.gamma_star)) + \
            np.sum((1-self.gamma_star) * np.log((1-self.gamma_star)))

        ELBO_loss = term1 + term2 - term3 + term4 - term5
        return ELBO_loss

    def fit(self, total_epochs=500):

        result = []
        for epoch in range(total_epochs):
            self.E_step()
            self.M_step()
            r_squared_test = self._evaluator.calculate_r_squared(
                self.mu_star, self._marginal_test, self._test_ld)
            result.append(r_squared_test)
            print("Epoch: {}, test R-squared: {}".format(epoch, r_squared_test))

        return result

    def predict(self, X):
        y_pred = X @ (self.gamma_star * self.mu_star)
        return y_pred

    def is_eligible(self, prior_type) -> bool:
        return prior_type == PRIOR_TYPE_SPIKE_AND_SLAB
