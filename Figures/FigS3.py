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
from matplotlib.patches import Patch


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



def suppfig3a(workdir, aligners, datasets):
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}
    matched_pcrt = []

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in ASMCALLERS:
            for aligner in ['minimap2', 'lra']:
                df_matched = pd.read_csv(f'{workdir}/{plat}/{dataset}_assembler_repro/ctg_assembler_repro/{aligner}-{caller}.assembler-concordant.tsv', sep='\t')
                df_unique = pd.read_csv(f'{workdir}/{plat}/{dataset}_assembler_repro/ctg_assembler_repro/{aligner}-{caller}.assembler-unique.tsv', sep='\t')
                assm = plat_assemblers[plat][0] + '-' + plat_assemblers[plat][1]
                matched_pcrt.append((caller, f'{assm} ({plat})', plat, 100 * len(df_matched) / (len(df_unique) + len(df_matched))))

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller in CALLERS:
            for i in range(len(aligners)):
                for j in range(len(aligners)):
                    if j <= i:
                        continue
                    df_matched = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/pw_aligner_repro/{caller}/{caller}.{aligners[i]}-{aligners[j]}.concordant.tsv', sep='\t')
                    df_unique = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/pw_aligner_repro/{caller}/{caller}.{aligners[i]}-{aligners[j]}.unique.tsv', sep='\t')
                    matched_pcrt.append((caller, f'{aligners[i]}-{aligners[j]}', plat, 100 * len(df_matched) / (len(df_unique) + len(df_matched))))

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_matched_pcrt = pd.DataFrame(matched_pcrt, columns=['caller', 'aligner', 'plat', 'pcrt'])

    my_pal = {aligner: '#fc8d62' if 'flye' in aligner else '#8da0cb' for aligner in df_matched_pcrt['aligner'].unique()}
    sns.violinplot(data=df_matched_pcrt, x='pcrt', y='aligner', palette=my_pal, scale='width', ax=ax)

    ax.set_xlabel('Percent of concordant SVs', fontsize=13)
    ax.set_xlim([25, 100])
    ax.set_xticks(np.linspace(25, 100, 4))
    ax.set_xticklabels([int(ele) for ele in np.linspace(25, 100, 4)], fontsize=12)
    ax.set_ylabel('')

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()

    matched_pcrt = []
    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in ASMCALLERS:
            for assembler in plat_assemblers[plat]:
                df_matched = pd.read_csv(
                    f'{workdir}/{plat}/{dataset}_aligner_repro/ctg_aligner_repro/{assembler}-{caller}.aligner-concordant.tsv',
                    sep='\t')
                df_unique = pd.read_csv(
                    f'{workdir}/{plat}/{dataset}_aligner_repro/ctg_aligner_repro/{assembler}-{caller}.aligner-unique.tsv',
                    sep='\t')
                matched_pcrt.append(
                    (caller, 'ctg-minimap2-lra', 100 * len(df_matched) / (len(df_matched) + len(df_unique))))

                if plat == 'ONT':
                    print(caller, assembler, 100 * len(df_matched) / (len(df_matched) + len(df_unique)))

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in CALLERS:
            for i in range(len(aligners)):
                for j in range(len(aligners)):
                    if j <= i:
                        continue
                    df_matched = pd.read_csv(
                        f'{workdir}/{plat}/{dataset}_aligner_repro/pw_aligner_repro/{caller}/{caller}.{aligners[i]}-{aligners[j]}.concordant.tsv',
                        sep='\t')
                    df_unique = pd.read_csv(
                        f'{workdir}/{plat}/{dataset}_aligner_repro/pw_aligner_repro/{caller}/{caller}.{aligners[i]}-{aligners[j]}.unique.tsv',
                        sep='\t')
                    matched_pcrt.append((caller, f'{aligners[i]}-{aligners[j]}',
                                         100 * len(df_matched) / (len(df_unique) + len(df_matched))))

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_matched_pcrt = pd.DataFrame(matched_pcrt, columns=['caller', 'aligner', 'pcrt'])

    my_pal = {aligner: '#fc8d62' if 'ctg' in aligner else '#8da0cb' for aligner in
              df_matched_pcrt['aligner'].unique()}
    sns.violinplot(data=df_matched_pcrt, x='pcrt', y='aligner', palette=my_pal, scale='width', ax=ax)

    ax.set_xlabel('Percent of concordant SVs', fontsize=13)
    ax.set_xlim([25, 100])
    ax.set_xticks(np.linspace(25, 100, 4))
    ax.set_xticklabels([int(ele) for ele in np.linspace(25, 100, 4)], fontsize=12)
    ax.set_ylabel('')

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    plt.show()


def suppfig3b(workdir, datasets):
    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'],
                     'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}
    unique_pcrt = []

    for caller in CALLERS:
        wgs_supp4_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}
        extd_supp4_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for dataset_idx, dataset in enumerate(datasets):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            this_dataset_index = datasets_dict[plat].index(dataset)
            unique_info = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.exTD.tsv', header=[0], sep='\t')
            matched_info = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-concordant.exTD.tsv', header=[0], sep='\t')

            merged_total = len(unique_info) + len(matched_info)

            supp_dict = {}
            for idx, row in matched_info.iterrows():
                supp = int(row['SUPP'])
                if supp in supp_dict:
                    supp_dict[supp] += 1
                else:
                    supp_dict[supp] = 1

            unique_pcrt.append((len(unique_info) / merged_total * 100, caller, plat, 'ExTD'))
            extd_supp4_pcrt[plat][this_dataset_index] += supp_dict[4] / merged_total * 100

            wgs_merged_vcf = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(wgs_merged_vcf)

            unique_pcrt.append((wgs_supp_dict[1] / wgs_merged_total * 100, caller, plat, 'WGS'))
            wgs_supp4_pcrt[plat][this_dataset_index] += wgs_supp_dict[4] / wgs_merged_total * 100


    df_unique = pd.DataFrame(unique_pcrt, columns=['pcrt', 'caller', 'plat', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    fig, axes = plt.subplots(1, 1, figsize=(3, 4))

    sns.barplot(data=df_unique, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                 palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[0])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of aligner unique SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['ExTD', 'WGS'], fontsize=13)

        ax.legend(handles=plat_legends)
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

        ax.grid(axis='y', ls='--', color='grey', lw=1.5)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()

    plt.show()

def suppfig3c(workdir, datasets, aligners):

    svtypes = ['INS', 'DEL', 'INV', 'DUP']
    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        fig, ax = plt.subplots(1, 1, figsize=(5, 4))

        svtypes_pcrt = []
        this_total = {aligner: 0 for aligner in aligners}
        svtypes_count = {aligner: {svtype: 0 for svtype in svtypes} for aligner in aligners}

        for i in range(len(aligners)):
            for j in range(len(aligners)):
                if j <= i:
                    continue
                for caller in CALLERS:
                    unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/pw_aligner_repro/{caller}/{caller}.{aligners[i]}-{aligners[j]}.unique.tsv'
                    df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])

                    for idx, row in df_unique.iterrows():
                        svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['ALIGNER'], \
                                                            row['REGION_TYPE']

                        this_total[aligner] += 1
                        if svtype in svtypes:
                            svtypes_count[aligner][svtype] += 1

        for aligner, region_counts in svtypes_count.items():
            for svtype, counts in region_counts.items():
                svtypes_pcrt.append((aligner, dataset, svtype, counts / this_total[aligner] * 100))

        df_svtypes_pcrt = pd.DataFrame(svtypes_pcrt, columns=['aligner', 'dataset', 'svtype', 'pcrt'])

        sns.barplot(data=df_svtypes_pcrt, x='svtype', y='pcrt', hue='aligner',
                    palette=['#1b9e77', '#d95f02', '#7570b3', '#e7298a'], ax=ax)

        ax.set_ylabel('% of read aligner specific SV types', fontsize=13)
        ax.set_ylim([0, 80])
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([int(ele) for ele in np.linspace(0, 80, 5)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['INS', 'DEL', 'INV', 'DUP'], fontsize=13)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        fig.tight_layout()

    plt.show()