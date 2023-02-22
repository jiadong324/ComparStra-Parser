#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/2/21

'''

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


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

def suppfig2a(workdir, aligners):

    plat_assembler = {'HiFi':['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plats = ['HiFi', 'ONT']

    all_supp_pcrt = {'HiFi': [], 'ONT': []}

    for fig_idx, plat in enumerate(plats):

        for caller in CALLERS:
            for aligner in aligners:
                merged_vcf = f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                for supp, count in supp_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, aligner, caller, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                merged_vcf = f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)

                for supp, count in supp_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, assembler, caller, 'Assembly'))

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    for col_idx, plat in enumerate(plats):
        this_ax = axes[col_idx]
        df_all_pcrt = pd.DataFrame(all_supp_pcrt[plat], columns=['pcrt', 'supp', 'aa', 'caller', 'stra'])
        sns.barplot(data=df_all_pcrt, x='supp', y='pcrt', hue='stra', hue_order=stra_order, palette=stra_color, capsize=.15, ax=this_ax)
        this_ax.set_title(plat, fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of WGS-SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.legend(handles=stra_legends)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.set_xlabel('# of datasets', fontsize=13)

    fig.tight_layout()

    all_supp_pcrt = {'HiFi': [], 'ONT': []}

    for fig_idx, plat in enumerate(plats):

        for caller in CALLERS:
            for aligner in aligners:
                matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv', header=[0], sep='\t')
                unique_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}.unique.exTD.tsv', header=[0], sep='\t')

                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1

                for supp, count in supps_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, aligner, caller, 'Read'))
                all_supp_pcrt[plat].append((len(unique_info) / merged_total * 100, 1, aligner, caller, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.{plat}-concordant.info.exTD.tsv',header=[0], sep='\t')
                unique_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.{plat}.unique.exTD.tsv', header=[0],sep='\t')
                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1

                for supp, count in supps_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, assembler, caller, 'Assembly'))
                all_supp_pcrt[plat].append((len(unique_info) / merged_total * 100, 1, assembler, caller, 'Assembly'))

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    for col_idx, plat in enumerate(plats):
        this_ax = axes[col_idx]
        df_all_pcrt = pd.DataFrame(all_supp_pcrt[plat], columns=['pcrt', 'supp', 'aa', 'caller', 'stra'])
        sns.barplot(data=df_all_pcrt, x='supp', y='pcrt', hue='stra', hue_order=stra_order,
                    palette=stra_color, capsize=.15, ax=this_ax)

        this_ax.set_title(plat, fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of ExTD-SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.legend(handles=stra_legends)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.set_xlabel('# of datasets', fontsize=13)

    fig.tight_layout()
    plt.show()

def suppfig2b(workdir, aligners):
    plat_assembler = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plats = ['HiFi', 'ONT']

    wgs_repro_pcrt = []
    extd_repro_pcrt = []

    for fig_idx, plat in enumerate(plats):

        for caller in CALLERS:
            for aligner in aligners:
                merged_vcf = f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                wgs_repro_pcrt.append((supp_dict[3] / merged_total * 100, aligner, TOOLMAP[caller], plat, 'Read'))

                matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv', header=[0], sep='\t')
                unique_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}.unique.exTD.tsv', header=[0], sep='\t')

                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1
                extd_repro_pcrt.append((supps_dict[3] / merged_total * 100, aligner, TOOLMAP[caller], plat, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                merged_vcf = f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                wgs_repro_pcrt.append((supp_dict[3] / merged_total * 100, assembler, TOOLMAP[caller], plat, 'Assembly'))

                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.{plat}-concordant.info.exTD.tsv',header=[0], sep='\t')
                unique_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.{plat}.unique.exTD.tsv', header=[0],sep='\t')
                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1
                extd_repro_pcrt.append((supps_dict[3] / merged_total * 100, assembler, TOOLMAP[caller], plat, 'Assembly'))

    df_wgs_repro = pd.DataFrame(wgs_repro_pcrt, columns=['pcrt', 'aa', 'caller', 'plat', 'stra'])
    df_extd_repro = pd.DataFrame(extd_repro_pcrt, columns=['pcrt', 'aa', 'caller', 'plat', 'stra'])

    read_wgs_repro, read_extd_repro = df_wgs_repro[df_wgs_repro['stra'] == 'Read'], df_extd_repro[df_extd_repro['stra'] == 'Read']
    assm_wgs_repro, assm_extd_repro = df_wgs_repro[df_wgs_repro['stra'] == 'Assembly'], df_extd_repro[df_extd_repro['stra'] == 'Assembly']

    assm_legends = [Line2D([0], [0], label='SVIM-asm', color=TOOLCOLORS['svimasm'], lw=2),
                   Line2D([0], [0], label='PAV', color=TOOLCOLORS['pav'], lw=2),
                   Line2D([0], [0], label='WGS', color='black', lw=2),
                   Line2D([0], [0], label='ExTD', color='black', ls='--', lw=2)]


    for idx, plat in enumerate(plats):
        fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 4))
        assm_plat_wgs, assm_plat_extd = assm_wgs_repro[assm_wgs_repro['plat'] == plat], assm_extd_repro[assm_extd_repro['plat'] == plat]

        sns.lineplot(data=assm_plat_wgs, x='aa', y='pcrt', hue='caller', hue_order=['PAV', 'SVIM-asm'], palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']],
                     lw=2, marker='o', markersize=8, ax=axes[0])

        sns.lineplot(data=assm_plat_extd, x='aa', y='pcrt', hue='caller', hue_order=['PAV', 'SVIM-asm'], palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']],
                     lw=2, marker='o', ls='--', markersize=8, ax=axes[1])

        for i, ax in enumerate(axes):

            if i == 0:
                ax.set_ylabel('% of datasets concordant SVs', fontsize=13)
            else:
                ax.set_ylabel('')
            ax.set_xlabel('')

            ax.set_xticks([0, 1])
            if plat == 'HiFi':
                ax.set_xticklabels(['hifiasm', 'flye'], fontsize=13, rotation=60)
            else:
                ax.set_xticklabels(['flye', 'shasta'], fontsize=13, rotation=60)
            ax.set_ylim(20, 100)
            ax.set_yticks(np.linspace(20, 100, 5))
            ax.set_yticklabels([int(val) for val in np.linspace(20, 100, 5)], fontsize=12)

            if i == 0:
                ax.legend(handles=assm_legends)
            else:
                ax.legend('', frameon=False)

            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        fig.tight_layout()

    for plat in plats:
        fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))

        sns.lineplot(data=read_wgs_repro[read_wgs_repro['plat'] == plat], x='aa', y='pcrt', hue='caller',
                    hue_order=[TOOLMAP[ele] for ele in CALLERS], palette=[TOOLCOLORS[ele] for ele in CALLERS],
                    lw=2, marker='o', markersize=8, ax=axes[0])

        sns.lineplot(data=read_extd_repro[read_extd_repro['plat'] == plat], x='aa', y='pcrt', hue='caller',
                    hue_order=[TOOLMAP[ele] for ele in CALLERS], palette=[TOOLCOLORS[ele] for ele in CALLERS],
                    lw=2, marker='o', ls='--', markersize=8, ax=axes[1])

        for idx, ax in enumerate(axes):
            if idx == 0:
                ax.set_ylabel('% of datasets concordant SVs', fontsize=13)
            else:
                ax.set_ylabel('')
            ax.set_xlabel('')

            ax.legend('', frameon=False)

            ax.set_xticks(np.arange(len(aligners)))
            ax.set_xticklabels(aligners, rotation=60, fontsize=13)

            ax.set_ylim(20, 100)
            ax.set_yticks(np.linspace(20, 100, 5))
            ax.set_yticklabels([int(val) for val in np.linspace(20, 100, 5)], fontsize=12)

            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        fig.tight_layout()

    plt.show()