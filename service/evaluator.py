import numpy as np


class Evaluator:

    def calculate_r_squared(self, mu_beta_s: np.ndarray,
                            marginal: np.ndarray,
                            ld: np.ndarray) -> float:

        result: float = None

        if mu_beta_s is not None and marginal is not None and ld is not None:

            ss_res = 1 - 2 * marginal.T @ mu_beta_s + \
                mu_beta_s.T @ ld @ mu_beta_s
            ss_total = 1

            result = 1 - ss_res / ss_total

        return result
