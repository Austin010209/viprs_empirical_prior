import argparse

from service.load_output_service import LoadOutputService
from service.plot_service import PlotService
from utils.constants import *

import numpy as np
import statistics
import os

"""
Main
"""

def fetch_data(chr):
    if not os.path.isdir(OUTPUT_DIR):
        raise FileNotFoundError(
            "Please create 'output' folder and train model first!")

    output_files = os.listdir(OUTPUT_DIR)
    if len(output_files) == 0:
        raise FileNotFoundError(
            "Please train model first!")

    print("Analyse chromosome {} output result".format(chr))

    # Load data
    load_output_service = LoadOutputService()
    result = load_output_service.load(output_files)

    if result is None or result.get(str(chr)) is None:
        raise FileNotFoundError(
            "Output result for chr: {} is not found!".format(chr))

    chr_output_result = result[str(chr)]

    r_squared_hist_empirical = chr_output_result.get("True")
    r_squared_hist_sas = chr_output_result.get("False")

    assert r_squared_hist_empirical is not None and len(
        r_squared_hist_empirical) > 0
    assert r_squared_hist_sas is not None and len(r_squared_hist_sas) > 0

    best_empirical = [np.amax(x) for x in r_squared_hist_empirical]
    best_sas = [np.amax(x) for x in r_squared_hist_sas]

    best_empirical_mean = statistics.mean(best_empirical)
    best_empirical_std = statistics.stdev(best_empirical)

    best_sas_mean = statistics.mean(best_sas)
    best_sas_std = statistics.stdev(best_sas)

    print("Emperical best r-squared: {:.3f}±{:.3f}".format(best_empirical_mean, best_empirical_std))
    print("Spike and slab best r-squared: {:.3f}±{:.3f}".format(best_sas_mean, best_sas_std))
    return r_squared_hist_empirical, r_squared_hist_sas, best_empirical, best_sas


def main(chr: int):
    r_squared_hist_empirical_21, r_squared_hist_sas_21, best_empirical_21, best_sas_21 = fetch_data(21)
    r_squared_hist_empirical_22, r_squared_hist_sas_22, best_empirical_22, best_sas_22 = fetch_data(22)

    # Plot
    plot_service = PlotService()
    plot_service.plot_r_squared_hist_all(TOTAL_FOLDS,
                                     TOTAL_EPOCHS,
                                     r_squared_hist_empirical_21, 
                                     r_squared_hist_sas_21,
                                     r_squared_hist_empirical_22, 
                                     r_squared_hist_sas_22,
                                     title="Test_R_squared_as_Training_Progresses_on_chromosome_21_and_22"
                                     )

    plot_service.plot_best_r_squared(best_empirical_21, 
                                     best_sas_21, 
                                     best_empirical_22, 
                                     best_sas_22,
                                     title="Test_R_squared_Comparison_between_Methods_on_Chromosome_21_and_22")


if __name__ == "__main__":

    # Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--chr', type=int,
                        choices=[21, 22],
                        default=22)
    args = parser.parse_args()

    main(args.chr)
