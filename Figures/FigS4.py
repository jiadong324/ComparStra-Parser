#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/2/22

'''

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math

from Helpers.Functions import *
from Helpers.Constant import *

sns.set_theme(style="ticks", font="Arial", font_scale=1.0)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2

def suppfig4a(workdir, dataset, region):

    for i, caller in enumerate(CALLERS):
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        wgs_type_dict = {}

        df_wgs = pd.read_csv(f'{workdir}/HiFi/read_callers_merged/type_repro/{dataset}_{caller}-concordant.info.{region}.tsv', sep='\t')
        for idx, row in df_wgs.iterrows():
            svtype = row['TYPE']
            if svtype in wgs_type_dict:
                wgs_type_dict[svtype] += 1
            else:
                wgs_type_dict[svtype] = 1

        wgs_log_num, wgs_pcrt = [], []
        xticks = []
        incorrect = 0
        for key, val in wgs_type_dict.items():
            if ',' in key:
                incorrect += val
            else:
                wgs_log_num.append(math.log10(val))
                wgs_pcrt.append(val / len(df_wgs))
                xticks.append(key)

        xticks.append('Disc')
        wgs_log_num.append(math.log10(incorrect))
        wgs_pcrt.append(incorrect / len(df_wgs))

        barwidth = 0.5

        bars = ax.bar(np.arange(len(xticks)), wgs_log_num, width=barwidth, color='#1b9e77')
        for j, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height + 0.01, f'{round(wgs_pcrt[j] * 100, 2)}%', ha='center', va='bottom')

        bars[-1].set_color('#e7298a')
        ax.set_xticks(np.arange(len(xticks)))
        ax.set_xticklabels(xticks, fontsize=13, rotation=60)
        ax.set_yticks([math.log10(1), math.log10(10), math.log10(100), math.log10(1000), math.log10(10000)])
        ax.set_yticklabels(['0', '10', '100', r'$10^3$', r'$10^4$'], fontsize=12)

        ax.set_ylabel(f'# of SVs from {caller}', fontsize=13)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        fig.tight_layout()
    plt.show()


def suppfig4b(workdir, dataset, aligners):

    for i, aligner in enumerate(aligners):
        fig, axes = plt.subplots(1, 2, figsize=(7, 4))
        wgs_type_dict = {}

        df_wgs = pd.read_csv(f'{workdir}/HiFi/read_callers_merged/type_repro/{dataset}_{aligner}-concordant.info.WGS.tsv', sep='\t')
        for idx, row in df_wgs.iterrows():
            svtype = row['TYPE']
            if svtype in wgs_type_dict:
                wgs_type_dict[svtype] += 1
            else:
                wgs_type_dict[svtype] = 1

        wgs_log_num, wgs_pcrt = [], []
        xticks = []
        incorrect = 0
        incorrect_dict = {'Rare': 0}
        for key, val in wgs_type_dict.items():
            if ',' in key:
                incorrect += val
                if val > 50:
                    incorrect_dict[key] = val
                else:
                    incorrect_dict['Rare'] += val
            else:
                wgs_log_num.append(math.log10(val))
                wgs_pcrt.append(val / len(df_wgs))
                xticks.append(key)

        xticks.append('Disc')
        wgs_log_num.append(math.log10(incorrect))
        wgs_pcrt.append(incorrect / len(df_wgs))

        barwidth = 0.5

        bars = axes[0].bar(np.arange(len(xticks)), wgs_log_num, width=barwidth, color = '#1b9e77')
        for j, bar in enumerate(bars):
            height = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width() / 2., height + 0.01, f'{round(wgs_pcrt[j] * 100, 2)}%', ha='center', va='bottom')

        bars[-1].set_color('#e7298a')
        axes[0].set_xticks(np.arange(len(xticks)))
        axes[0].set_xticklabels(xticks, fontsize=13)
        axes[0].set_yticks([math.log10(1), math.log10(10), math.log10(100), math.log10(1000), math.log10(10000)])
        axes[0].set_yticklabels(['0', '10', '100', r'$10^3$', r'$10^4$'], fontsize=12)


        axes[0].set_ylabel('NO. of SVs', fontsize=13)

        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)


        sorted_incorrect = sorted(incorrect_dict.items(), key=lambda x:x[1])
        bars = axes[1].barh(np.arange(len(sorted_incorrect)), [math.log10(ele[1]) for ele in sorted_incorrect], height=barwidth, color='#377eb8')
        bars[-1].set_color('#d95f02')
        axes[1].set_yticks(np.arange(len(sorted_incorrect)))
        axes[1].set_yticklabels([ele[0] for ele in sorted_incorrect])
        axes[1].set_xticks([math.log10(1), math.log10(10), math.log10(100), math.log10(1000)])
        axes[1].set_xticklabels(['0', '10', '100', r'$10^3$'], fontsize=12)
        axes[1].set_xlabel('NO. of SVs', fontsize=13)
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['right'].set_visible(False)

        fig.tight_layout()

    plt.show()