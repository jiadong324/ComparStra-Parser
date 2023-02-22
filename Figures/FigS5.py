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


def suppfig5a(workdir, datasets):
    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'],
                     'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}

    fig, axes = plt.subplots(1, 1, figsize=(3, 4))
    unique_pcrt = []

    for caller in ASMCALLERS:
        wgs_supp2_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}
        extd_supp2_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for dataset_idx, dataset in enumerate(datasets):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            unique_info = pd.read_csv(f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.assembler-unique.exTD.tsv', header=[0], sep='\t')
            matched_info = pd.read_csv(f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.assembler-concordant.exTD.tsv', header=[0], sep='\t')

            merged_total = len(unique_info) + len(matched_info)
            supp_dict = {}

            for idx, row in matched_info.iterrows():
                supp = int(row['SUPP'])
                if supp in supp_dict:
                    supp_dict[supp] += 1
                else:
                    supp_dict[supp] = 1

            this_dataset_index = datasets_dict[plat].index(dataset)
            unique_pcrt.append((len(unique_info) / merged_total * 100, plat, caller, 'ExTD'))

            extd_supp2_pcrt[plat][this_dataset_index] += supp_dict[2] / merged_total * 100

            merged_vcf = f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(merged_vcf)

            unique_pcrt.append((wgs_supp_dict[1] / wgs_merged_total * 100, plat, caller, 'WGS'))
            wgs_supp2_pcrt[plat][this_dataset_index] += wgs_supp_dict[2] / wgs_merged_total * 100

    df_unique = pd.DataFrame(unique_pcrt, columns=['pcrt', 'plat', 'caller', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    sns.barplot(data=df_unique, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[0])


    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of assembler unique SVs', fontsize=13)
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

def suppfig5b(workdir, datasets):
    regions = ['Repeat Masked', 'Segment Dup', 'Simple Region', 'Tandem Repeats']
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'],
                       'ONT': ['flye', 'shasta']}

    svtypes_info = []

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(ASMCALLERS):
            unique_info_out = f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.assembler-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            svtypes_dict = {assm: {'INS': 0, 'DEL': 0} for assm in plat_assemblers[plat]}
            this_total = {assm: 0 for assm in plat_assemblers[plat]}
            assembler_unique_region_dict = {assembler: [0 for i in range(len(regions))] for assembler in
                                            plat_assemblers[plat]}
            assembler_unique_dict = {assembler: 0 for assembler in plat_assemblers[plat]}

            for idx, row in df_unique.iterrows():
                assembler, svtype, svlen, sv_region = row['ASSEMBLER'], row['TYPE_MATCH'], abs(int(row['SVLEN'])), row[
                    'REGION_TYPE']
                if svtype in ['INS', 'DEL']:
                    svtypes_dict[assembler][svtype] += 1
                    this_total[assembler] += 1
                    region_idx = regions.index(sv_region)
                    assembler_unique_dict[assembler] += 1
                    assembler_unique_region_dict[assembler][region_idx] += 1

            for assm, svtype_count in svtypes_dict.items():
                for svtype, count in svtype_count.items():
                    svtypes_info.append((svtype, count / this_total[assm] * 100, assm, plat, dataset, caller))


    df_svtypes_info = pd.DataFrame(svtypes_info, columns=['svtype', 'pcrt', 'assembler', 'plat', 'dataset', 'caller'])

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    sns.boxplot(data=df_svtypes_info[df_svtypes_info['plat'] == 'HiFi'], x='svtype', y='pcrt', hue='assembler',
                hue_order=['flye', 'hifiasm'],
                palette=[ASSMBLERCOLOR['flye'], ASSMBLERCOLOR['hifiasm']], ax=axes[0])

    sns.boxplot(data=df_svtypes_info[df_svtypes_info['plat'] == 'ONT'], x='svtype', y='pcrt', hue='assembler',
                hue_order=['flye', 'shasta'], palette=[ASSMBLERCOLOR['flye'], ASSMBLERCOLOR['shasta']], ax=axes[1])

    for i, ax in enumerate(axes):
        ax.set_ylim(0, 90)
        ax.set_yticks(np.linspace(0, 90, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 90, 4)], fontsize=12)
        ax.legend(title='Assembler', loc='lower left')
        if i == 0:
            ax.set_ylabel('% of assembler unique SV', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)
        ax.set_xlabel('')

    fig.tight_layout()
    plt.show()

def suppfig5c(workdir, datasets):
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        aligner_unique_types = []

        for caller in ASMCALLERS:
            for assembler in plat_assemblers[plat]:
                df_uniques = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/ctg_aligner_repro/{assembler}-{caller}.aligner-unique.tsv', sep='\t')
                unique_types_dict = {'lra':{}, 'minimap2': {}}
                unique_total = {'lra': 0, 'minimap2': 0}
                for idx, row in df_uniques.iterrows():
                    svtype, aligner = row['TYPE_MATCH'], row['ALIGNER']
                    unique_total[aligner] += 1
                    if svtype in unique_types_dict[aligner]:
                        unique_types_dict[aligner][svtype] += 1
                    else:
                        unique_types_dict[aligner][svtype] = 1

                for aligner, count in unique_types_dict.items():
                    for key, val in count.items():
                        aligner_unique_types.append((aligner, key, val / unique_total[aligner] * 100, caller, assembler))


        df_aligner_unique_type = pd.DataFrame(aligner_unique_types, columns=['aligner', 'svtype', 'pcrt', 'caller', 'assembler'])

        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        sns.barplot(data=df_aligner_unique_type, x='svtype', y='pcrt', hue='aligner', hue_order=['lra', 'minimap2'],
                    palette=['#1b9e77', '#d95f02'], capsize=.2, ax=ax)

        ax.set_ylabel('% of contig aligner specific SV types', fontsize=13)
        ax.set_ylim([0, 80])
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([int(ele) for ele in np.linspace(0, 80, 5)], fontsize=12)


        ax.set_xlabel('')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['INS', 'DEL', 'INV'], fontsize=13)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

        fig.tight_layout()

        plt.show()
