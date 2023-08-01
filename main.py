import argparse

from model.viprs_factory import ViprsFactory
from loader.chr_data_loader import ChrDataLoader
from service.evaluator import Evaluator
from utils.constants import *

import numpy as np
import os

"""
Main
"""


def main(use_annot: bool, chr: int):

    prior_type = PRIOR_TYPE_EMPIRICAL if use_annot else PRIOR_TYPE_SPIKE_AND_SLAB

    data_loader = ChrDataLoader()
    evaluator = Evaluator()

    for fold in range(TOTAL_FOLDS):

        fold += 1
        print("Fold: {}".format(fold))

        train_data, test_data = data_loader.load(chr, fold)
        factory = ViprsFactory(train_data, test_data, evaluator,
                               learning_rate=LEARNING_RATE,
                               tau_beta=TAU_BETA_INIT,
                               tau_eps_sas=TAU_EPS_INIT_SAS,
                               tau_eps_empirical=TAU_EPS_INIT_EMPIRICAL,
                               tau_beta_s=TAU_BETA_S_INIT,
                               mu_beta_s=MU_BETA_S_INIT,
                               gamma_s=GAMMA_S_INIT,
                            #    var_per_snp=VAR_PER_SNP,
                               tau_beta_init_empirical=TAU_BETA_INIT_EMPIRICAL,
                               pi_sas=PI_INIT_SAS,
                               pi_empirical=PI_INIT_EMPIRICAL)

        viprs = factory.get_viprs(prior_type)
        test_r_squared = viprs.fit(TOTAL_EPOCHS)
        np.savetxt('{}/fold_{}_chr_{}_use_annot_{}.csv'.format(OUTPUT_DIR,
                   fold, chr, use_annot), test_r_squared)


if __name__ == "__main__":

    # Check directory
    if not os.path.isdir(DATA_DIR):
        raise FileNotFoundError(
            "Please create 'data' folder and download data first!")

    if not os.path.isdir(OUTPUT_DIR):
        raise FileNotFoundError(
            "Please create 'output' folder first!")

    # Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--use_annot', type=str,
                        choices=['true', 'false'],
                        default='true')
    parser.add_argument('--chr', type=int,
                        choices=[21, 22],
                        default=22)
    args = parser.parse_args()

    use_annot = args.use_annot.lower() == 'true'
    if use_annot:
        print('Use annotation (empirical prior VIPRS)')
    else:
        print('Use baseline (spike and slab prior VIPRS)')

    main(use_annot, args.chr)
