import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List

from utils.constants import *


class PlotService:

    def plot_r_squared_hist(self, total_folds: int,
                            total_epochs: int,
                            r_squared_hist_empirical: List[np.ndarray],
                            r_squared_hist_sas: List[np.ndarray],
                            title: str) -> None:

        history_list = []
        for fold in range(total_folds):

            r_squared_empirical = r_squared_hist_empirical[fold]
            r_squared_sas = r_squared_hist_sas[fold]

            epochs = list(range(total_epochs))
            trial = [str(fold) for _ in range(total_epochs)]

            data_preproc = pd.DataFrame({
                'Epochs': epochs,
                'Fold': trial,
                'Use annotations': r_squared_empirical,
                'Baseline': r_squared_sas})

            history = pd.melt(data_preproc, id_vars=['Epochs', 'Fold'])
            history_list.append(history)

        r_squared_df = pd.concat(history_list, ignore_index=True, sort=False)
        r_squared_df.columns = ["Epochs", "Fold", "Algorithm", "R-squared"]

        plt.clf()
        plt.grid(linewidth=0.5)
        sns.lineplot(x='Epochs', y='R-squared', hue='Algorithm',
                     data=r_squared_df).set(title=title)
        plt.savefig("{}/{}.png".format(OUTPUT_DIR, title))
        plt.close()
    

    def get_data_ready_1(self, total_folds, total_epochs, r_squared_hist_empirical, r_squared_hist_sas):
        history_list = []
        for fold in range(total_folds):

            r_squared_empirical = r_squared_hist_empirical[fold]
            r_squared_sas = r_squared_hist_sas[fold]

            epochs = list(range(total_epochs))
            trial = [str(fold) for _ in range(total_epochs)]

            data_preproc = pd.DataFrame({
                'Epochs': epochs,
                'Fold': trial,
                'Use annotations': r_squared_empirical,
                'Baseline': r_squared_sas})

            history = pd.melt(data_preproc, id_vars=['Epochs', 'Fold'])
            history_list.append(history)

        r_squared_df = pd.concat(history_list, ignore_index=True, sort=False)
        r_squared_df.columns = ["Epochs", "Fold", "Algorithm", "R-squared"]
        return r_squared_df
       
    

    def plot_r_squared_hist_all(self, total_folds: int,
                            total_epochs: int,
                            r_squared_hist_empirical_21, 
		                     r_squared_hist_sas_21,
		                     r_squared_hist_empirical_22, 
		                     r_squared_hist_sas_22,
                            title: str) -> None:
        r_squared_df_21 = self.get_data_ready_1(total_folds, total_epochs, r_squared_hist_empirical_21, r_squared_hist_sas_21)
        r_squared_df_22 = self.get_data_ready_1(total_folds, total_epochs, r_squared_hist_empirical_22, r_squared_hist_sas_22)
        

        plt.clf()
        fig, axes = plt.subplots(1, 2, sharey = True, figsize=(10, 6))
        ax0 = axes[0]
        ax0.grid(linewidth=0.5)
        ax0.set_title("Chromosome 21")
        sns.lineplot(ax=ax0, x='Epochs', y='R-squared', hue='Algorithm',
                     data=r_squared_df_21)#.set(title=title)
        ax1 = axes[1]
        ax1.grid(linewidth=0.5)
        ax1.set_title("Chromosome 22")
        sns.lineplot(ax=ax1, x='Epochs', y='R-squared', hue='Algorithm',
                     data=r_squared_df_22)#.set(title=title)
        fig.suptitle(title)
        plt.savefig("{}/final_figures/{}.png".format(OUTPUT_DIR, title))
        plt.close()


    def plot_best_r_squared(self, best_empirical_21: List[float],
                            best_sas_21: List[float],
                            best_empirical_22: List[float],
                            best_sas_22: List[float],
                            title: str) -> None:

        best_r_squared_21 = [best_empirical_21, best_sas_21]
        best_r_squared_22 = [best_empirical_22, best_sas_22]

        plt.clf()
        fig, axes = plt.subplots(1, 2, sharey = True, figsize=(10, 6))
        ax0 = axes[0]
        ax0.boxplot(best_r_squared_21,
                    labels=['Use annotations', 'Baseline'])
        ax0.set_xlabel('Algorithm')
        ax0.set_ylabel('R-squared')
        ax0.set_title("Chromosome 21")
        
        ax1 = axes[1]
        ax1.boxplot(best_r_squared_22,
                    labels=['Use annotations', 'Baseline'])
        ax1.set_xlabel('Algorithm')
        ax1.set_title("Chromosome 22")
        # ax1.set_ylabel('R-squared')

        fig.suptitle(title)
        plt.savefig("{}/final_figures/{}.png".format(OUTPUT_DIR, title))
        plt.close()
