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

def figure2a(workdir, datasets, aligners):

    ins_pcrt = []
    del_pcrt = []

    sv_count_list = []
    align_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in align_sv_count.iterrows():
        plat = 'HiFi'
        all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_num']), int(row['ins_num']), int(
            row['del_num']), row['aligner'], row['dataset'], row['caller']
        if 'ont' in dataset:
            plat = 'ONT'

        ins_pcrt.append((TOOLMAP[caller], aligner, PLATMAP[dataset], ins_num * 100 / all_sv_num, plat, 'Read-based'))
        # del_pcrt.append((TOOLMAP[caller], aligner, PLATMAP[dataset], del_num * 100 / all_sv_num, plat, 'Read-based'))
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], all_sv_num, plat, 'Read-based'))

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        svcount = ins_num + del_num + int(row['inv_num'])
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))
        # del_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], del_num * 100 / svcount, plat, 'Assembly-based'))

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))

    df_ins_pcrt = pd.DataFrame(ins_pcrt, columns=['caller', 'aa', 'dataset', 'pcrt', 'plat', 'stra'])
    df_sv_counts = pd.DataFrame(sv_count_list, columns=['caller', 'dataset', 'count', 'plat', 'stra'])

    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]
    stra_legend = [Patch(label='Assembly', color=STRACOLORS['Assembly']), Patch(label='Read', color=STRACOLORS['Read'])]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))
    sns.boxplot(data=df_sv_counts, y='dataset', x='count', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[0])
    sns.boxplot(data=df_ins_pcrt, y='dataset', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[1])
    for i, ax in enumerate(axes):

        ax.set_ylabel('')
        ax.set_yticks([0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb', 'ONT-9kb', 'ONT-19kb', 'ONT-30kb'], fontsize=13)
        if i == 0:
            ax.set_xlim(10000, 40000)
            ax.set_xticks(np.linspace(10000, 40000, 4))
            ax.set_xticklabels([int(val) for val in np.linspace(10, 40, 4)], fontsize=12)

            ax.set_xlabel('# of SVs (x$10^3$)', fontsize=13)
        else:
            ax.set_xlim(0, 100)
            ax.set_xticks(np.linspace(0, 100, 5))
            ax.set_xticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

            ax.set_xlabel('% of INS', fontsize=13)

        ax.legend(handles=stra_legend)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='x', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    plt.show()

def figure2b(workdir, aligners):

    plat_assembler = {'HiFi':['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plats = ['HiFi', 'ONT']

    supp_pcrt = []

    for fig_idx, plat in enumerate(plats):

        for caller in CALLERS:
            for aligner in aligners:
                merged_vcf = f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                supp_pcrt.append((supp_dict[3] / merged_total, aligner, TOOLMAP[caller], plat, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                merged_vcf = f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                supp_pcrt.append((supp_dict[3] / merged_total, assembler, TOOLMAP[caller], plat, 'Assembly'))

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]


    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'aa', 'caller', 'plat', 'stra'])

    sns.barplot(data=df_supp_pcrt, x='plat', y='pcrt', hue='stra',
                       hue_order=stra_order, palette=stra_color, capsize=0.15, ax=ax)

    ax.set_ylabel('% of concordant SVs', fontsize=13)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    ax.legend(handles=stra_legends)

    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels([f'{int(val * 100)}' for val in np.linspace(0, 1, 5)], fontsize=12)

    ax.set_xlabel('')
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(['HiFi', 'ONT'], fontsize=13)

    fig.tight_layout()
    plt.show()


def figure2c(workdir, aligners):

    shift_labels = ['0', '0,10', '>10']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']
    bpstd = []

    plat_assembler = {'hifi':['hifiasm', 'flye'], 'ont':['shasta', 'flye']}

    for k, plat in enumerate(['hifi', 'ont']):

        wgs_start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
        extd_start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}

        wgs_caller_total = [1 for k in range(len(callers))]
        extd_caller_total = [1 for k in range(len(callers))]


        for asm_caller in ASMCALLERS:
            caller_index = callers.index(asm_caller)
            for assembler in plat_assembler[plat]:

                wgs_matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{PLATMAP[plat]}/{asm_caller}.{assembler}.{PLATMAP[plat]}-concordant.info.tsv', sep='\t',header=[0])
                for idx, row in wgs_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        wgs_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            wgs_start_std_dict['0'][caller_index] += 1

                extd_matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{PLATMAP[plat]}/{asm_caller}.{assembler}.{PLATMAP[plat]}-concordant.info.exTD.tsv', sep='\t', header=[0])
                for idx, row in extd_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        extd_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            extd_start_std_dict['0'][caller_index] += 1


        for caller in CALLERS:
            caller_index = callers.index(caller)
            for aligner in aligners:
                wgs_matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}-concordant.info.tsv',sep='\t', header=[0])
                for idx, row in wgs_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        wgs_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            wgs_start_std_dict['0'][caller_index] += 1

                extd_matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv', sep='\t', header=[0])
                for idx, row in extd_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        extd_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            extd_start_std_dict['0'][caller_index] += 1

        extd_start_shift = [i / j * 100 for i, j in zip(extd_start_std_dict['0'], extd_caller_total)]
        for i, pcrt in enumerate(extd_start_shift):
            if i < len(CALLERS):
                bpstd.append((pcrt, PLATMAP[plat], callers[i], 'Read', 0, 'ExTD'))
            else:
                bpstd.append((pcrt, PLATMAP[plat], callers[i], 'Assembly', 0, 'ExTD'))

        wgs_start_shift = [i / j * 100 for i, j in zip(wgs_start_std_dict['0'], wgs_caller_total)]

        for i, pcrt in enumerate(wgs_start_shift):
            if i < len(CALLERS):
                bpstd.append((pcrt, PLATMAP[plat], callers[i], 'Read', 0, 'WGS'))
            else:
                bpstd.append((pcrt, PLATMAP[plat], callers[i], 'Assembly', 0, 'WGS'))



    df_bpstd = pd.DataFrame(bpstd, columns=['pcrt', 'plat', 'caller', 'stra', 'shift', 'region'])
    # df_extd_bpstd = pd.DataFrame(extd_bpstd, columns=['pcrt', 'plat', 'caller', 'shift', 'stra'])
    hifi_bpstd = df_bpstd[df_bpstd['plat'] == 'HiFi']
    ont_bpstd = df_bpstd[df_bpstd['plat'] == 'ONT']


    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 4))

    sns.boxplot(data=hifi_bpstd, x='region', y='pcrt', hue='stra', hue_order=['Read', 'Assembly'],
                palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], ax=axes[0])
    axes[0].set_title('HiFi', fontsize=13)

    sns.boxplot(data=ont_bpstd, x='region', y='pcrt', hue='stra', hue_order=['Read', 'Assembly'],
                palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], ax=axes[1])
    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of BSD-0 SVs', fontsize=13)
            ax.legend(title='')
        else:
            ax.set_ylabel('')
            ax.legend('', frameon=False)

        ax.set_xlabel('')
        # ax.set_xticks([0, 1])
        # ax.set_xticklabels(['WGS', 'ExTD'], fontsize=13)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()


    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    barwidth = 0.35
    hifi_xticks = np.arange(len(callers))
    ont_xticks = [r + barwidth + 0.06 for r in hifi_xticks]
    xticks = [r + (barwidth + 0.06) / 2 for r in hifi_xticks]

    rs = [hifi_xticks, ont_xticks]

    for plat_idx, plat in enumerate(['HiFi', 'ONT']):
        shift_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
        totals = [0 for k in range(len(callers))]

        for caller in CALLERS:
            caller_index = callers.index(caller)
            for row_idx, aligner in enumerate(aligners):
                matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}-concordant.info.tsv', sep='\t',header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region, rep_pcrt = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE'], float(row['PCRT'])

                    totals[caller_index] += 1

                    if start_std == 0:
                        shift_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        shift_dict['0,10'][caller_index] += 1
                    else:
                        shift_dict['>10'][caller_index] += 1

        assemblers = plat_assembler[plat]
        for caller in ASMCALLERS:
            caller_index = callers.index(caller)
            for assembler in assemblers:
                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{plat}/{caller}.{assembler}.{plat}-concordant.info.tsv', sep='\t', header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), \
                                                                   row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                    totals[caller_index] += 1

                    if start_std == 0:
                        shift_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        shift_dict['0,10'][caller_index] += 1
                    else:
                        shift_dict['>10'][caller_index] += 1

        for l, shift in enumerate(shift_labels):

            this_shift = [i / j * 100 for i, j in zip(shift_dict[shift], totals)]

            if l == 0:
                if plat == 'HiFi':
                    ax.bar(rs[plat_idx], this_shift, color=SHIFTCOLORDICT[shift_labels[l]], width=barwidth,
                                edgecolor='black', label=shift)
                else:
                    ax.bar(rs[plat_idx], this_shift, color=SHIFTCOLORDICT[shift_labels[l]], width=barwidth,
                                edgecolor='black', label=shift, hatch='\\\\')
            else:
                bottom = []
                for m in range(0, l):
                    bottom.append([i / j * 100 for i, j in zip(shift_dict[shift_labels[m]], totals)])
                bottom_sum = [sum(x) for x in zip(*bottom)]
                if plat == 'HiFi':
                    ax.bar(rs[plat_idx], this_shift, bottom=bottom_sum, color=SHIFTCOLORDICT[shift_labels[l]], width=barwidth, edgecolor='black', label=shift)
                else:
                    ax.bar(rs[plat_idx], this_shift, bottom=bottom_sum, color=SHIFTCOLORDICT[shift_labels[l]], width=barwidth, edgecolor='black', label=shift, hatch='\\\\')

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)

    ax.set_xlabel('')
    ax.set_xticks(xticks)
    ax.set_xticklabels([TOOLMAP[caller] for caller in callers], fontsize=13)

    ax.legend('', frameon=False)

    ax.set_ylabel('% of concordant SVs', fontsize=13)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

    fig.tight_layout()
    plt.show()

def figure2d(workdir, datasets):
    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'],
                     'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}
    concordants_pcrt = []

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

            concordants_pcrt.append((supp_dict[4] / merged_total * 100, caller, plat, 'ExTD'))
            extd_supp4_pcrt[plat][this_dataset_index] += supp_dict[4] / merged_total * 100

            wgs_merged_vcf = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(wgs_merged_vcf)

            concordants_pcrt.append((wgs_supp_dict[4] / wgs_merged_total * 100, caller, plat, 'WGS'))
            wgs_supp4_pcrt[plat][this_dataset_index] += wgs_supp_dict[4] / wgs_merged_total * 100


    df_concordants = pd.DataFrame(concordants_pcrt, columns=['pcrt', 'caller', 'plat', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    fig, axes = plt.subplots(1, 1, figsize=(3, 4))

    sns.barplot(data=df_concordants, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                 palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of aligner concordant SVs', fontsize=13)
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

    fig, axes = plt.subplots(1, 1, figsize=(3, 4))
    unique_pcrt = []
    concordants_pcrt = []

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
            concordants_pcrt.append((supp_dict[2] / merged_total * 100, plat, caller, 'ExTD'))

            extd_supp2_pcrt[plat][this_dataset_index] += supp_dict[2] / merged_total * 100

            merged_vcf = f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(merged_vcf)

            unique_pcrt.append((wgs_supp_dict[1] / wgs_merged_total * 100, plat, caller, 'WGS'))
            concordants_pcrt.append((wgs_supp_dict[2] / wgs_merged_total * 100, plat, caller, 'WGS'))
            wgs_supp2_pcrt[plat][this_dataset_index] += wgs_supp_dict[2] / wgs_merged_total * 100

    df_concordants = pd.DataFrame(concordants_pcrt, columns=['pcrt', 'plat', 'caller', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    sns.barplot(data=df_concordants, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of aligner concordant SVs', fontsize=13)
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

def aligner_assembler_wgs_extd_bpstd(workdir, datasets):
    leftbp_pcrt_list = []
    shift_labels = ['0', '0,10']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    for platform_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller_idx, caller in enumerate(CALLERS):
            wgs_shift_label = {shift: 0 for shift in shift_labels}
            df_matched = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-concordant.tsv', sep='\t', header=[0])
            this_total = 0

            for idx, row in df_matched.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 4:
                    this_total += 1
                    if start_std == 0:
                        wgs_shift_label['0'] += 1
                    elif start_std <= 10:
                        wgs_shift_label['0,10'] += 1

            df_extd_match_info = pd.read_csv(f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-concordant.exTD.tsv', sep='\t', header=[0])
            extd_total = 0
            extd_shift_label = {shift: 0 for shift in shift_labels}
            for idx, row in df_extd_match_info.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 4:
                    extd_total += 1
                    if start_std == 0:
                        extd_shift_label['0'] += 1
                    elif start_std <= 10:
                        extd_shift_label['0,10'] += 1

            for shift in shift_labels:
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller], wgs_shift_label[shift] / this_total * 100, shift, 'WGS'))
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller],extd_shift_label[shift] / extd_total * 100, shift, 'ExTD'))

        for caller in ASMCALLERS:
            wgs_shift_label = {shift: 0 for shift in shift_labels}
            df_matched = pd.read_csv(f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.assembler-concordant.tsv', sep='\t', header=[0])
            this_total = 0

            for idx, row in df_matched.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 2:
                    this_total += 1
                    if start_std == 0:
                        wgs_shift_label['0'] += 1
                    elif start_std <= 10:
                        wgs_shift_label['0,10'] += 1


            df_extd_matched = pd.read_csv(f'{workdir}/{plat}/{dataset}_assembler_repro/{caller}.assembler-concordant.exTD.tsv', sep='\t', header=[0])
            extd_total = 0
            extd_shift_label = {shift: 0 for shift in shift_labels}
            for idx, row in df_extd_matched.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 2:
                    extd_total += 1
                    if start_std == 0:
                        extd_shift_label['0'] += 1
                    elif start_std <= 10:
                        extd_shift_label['0,10'] += 1

            for shift in shift_labels:
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller], wgs_shift_label[shift] / this_total * 100, shift, 'WGS'))
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller], extd_shift_label[shift] / extd_total * 100, shift, 'ExTD'))

    df_bpstd = pd.DataFrame(leftbp_pcrt_list, columns=['plat', 'dataset', 'caller', 'pcrt', 'shift', 'region'])
    df_bpstd_subset1 = df_bpstd[df_bpstd['shift'] == '0']

    caller_legends = [Line2D([0], [0], color='white', mfc=TOOLCOLORS[ele], marker='o', markersize=8, label=TOOLMAP[ele]) for ele in callers]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
    sns.stripplot(data=df_bpstd_subset1[df_bpstd_subset1['region'] == 'WGS'], x='dataset', y='pcrt', hue='caller', marker='o', lw=2, size=8,
                  hue_order=[TOOLMAP[ele] for ele in callers], palette=[TOOLCOLORS[ele] for ele in callers], ax=axes[0])
    sns.stripplot(data=df_bpstd_subset1[df_bpstd_subset1['region'] == 'ExTD'], x='dataset', y='pcrt', hue='caller', marker='o', lw=2, size=8,
                  hue_order=[TOOLMAP[ele] for ele in callers], palette=[TOOLCOLORS[ele] for ele in callers], ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of BSD-0 SVs', fontsize=13)
            ax.legend(handles=caller_legends)
        else:
            ax.set_ylabel('')
            ax.legend('', frameon=False)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks(np.arange(len(datasets)))
        ax.set_xticklabels([PLATMAP[ele] for ele in datasets], fontsize=13, rotation=60)


        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    plt.show()

def figure2f(workdir, datasets, aligners):
    aligner_unique_counts = []
    aligner_unique_info = []
    aligner_uniques = {aligner: [] for aligner in aligners}

    unique_size_in_range = {aligner: 0 for aligner in aligners}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for j, caller in enumerate(CALLERS):

            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.exTD.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            aligner_unique_dict = {aligner: 0 for aligner in aligners}
            aligner_unique_region_dict = {aligner: [0 for i in range(len(GENOMICREGIONS))] for aligner in aligners}

            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row[
                    'REGION_TYPE']
                region_idx = GENOMICREGIONS.index(sv_region)
                aligner_unique_dict[aligner] += 1
                aligner_unique_region_dict[aligner][region_idx] += 1
                if svlen >= 50:
                    aligner_unique_info.append((svlen, aligner, sv_region, svtype, plat))

                if svlen >= 100 and svlen <= 1000:
                    unique_size_in_range[aligner] += 1

            for aligner, count in aligner_unique_dict.items():
                aligner_unique_counts.append((count, plat, PLATMAP[dataset], caller, aligner))
                aligner_uniques[aligner].append(count)

    for aligner, counts in aligner_uniques.items():
        print(f'{aligner} median: {np.median(counts)} avg: {np.mean(counts)}')

    print(unique_size_in_range)

    df_aligner_unique = pd.DataFrame(aligner_unique_info, columns=['svlen', 'aligner', 'region', 'svtype', 'plat'])
    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(6, 4))

    aligner_order = ['winnowmap', 'lra', 'ngmlr', 'minimap2']
    aligner_legends = [Patch(facecolor=ALIGNERCOLOR[aligner], label=aligner) for aligner in aligner_order]

    sns.histplot(data=df_aligner_unique[df_aligner_unique['plat'] == 'HiFi'], x='svlen', hue='aligner', bins=50,
                 alpha=1,hue_order=['minimap2', 'winnowmap', 'lra', 'ngmlr'], palette=[ALIGNERCOLOR[ele] for ele in ['minimap2', 'winnowmap', 'lra', 'ngmlr']], log_scale=True,
                 ax=axes[0])

    axes[0].set_ylabel('HiFi (x$10^3$)', fontsize=13)

    sns.histplot(data=df_aligner_unique[df_aligner_unique['plat'] == 'ONT'], x='svlen', hue='aligner', bins=50, alpha=1,
                 hue_order=['winnowmap', 'minimap2', 'ngmlr', 'lra'], palette=[ALIGNERCOLOR[ele] for ele in ['winnowmap', 'minimap2', 'ngmlr', 'lra']], log_scale=True, ax=axes[1])

    axes[1].set_ylabel('ONT (x$10^3$)', fontsize=13)

    for ax in axes:
        ax.legend(handles=aligner_legends)
        ax.tick_params(axis='both', which='both', length=0)

        ax.set_xlabel('SV length(bp)', fontsize=12)
        ax.set_xlim(50, 100000)
        ax.set_xticks([50, 100, 300, 1000, 6000, 100000])
        ax.set_xticklabels(['50', '100', '300', '1,000', '6,000', '100,000'], fontsize=12, rotation=60)

        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_ylabel('Number of SVs', fontsize=13)

        ax.set_xlabel('')
        # ax.set_ylim(0, 1000)
        ax.set_yticks(np.linspace(0, 1000, 3))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 1000, 3)], fontsize=12)

    fig.tight_layout()
    plt.show()

def figure2g(workdir, datasets, aligners):

    svtypes = ['INS', 'DUP', 'DEL']

    svtypes_pcrt = {aligner: {region: {svtype: [] for svtype in svtypes} for region in GENOMICREGIONS} for aligner in aligners}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(CALLERS):
            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.exTD.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])

            this_total = {aligner: 0 for aligner in aligners}
            svtypes_count = {aligner: {region: [0, 0, 0] for region in GENOMICREGIONS} for aligner in aligners}
            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']

                if svlen >= 100 and svlen <= 1000:
                    this_total[aligner] += 1
                    if svtype in svtypes:
                        svtype_idx = svtypes.index(svtype)
                        svtypes_count[aligner][sv_region][svtype_idx] += 1

            for aligner, region_counts in svtypes_count.items():
                for region, counts in region_counts.items():
                    for i, val in enumerate(counts):
                        svtype = svtypes[i]
                        svtypes_pcrt[aligner][region][svtype].append(val / this_total[aligner] * 100)

    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(7, 4))
    xticks = np.arange(len(svtypes))
    bar_width = 0.6

    for col_idx, aligner in enumerate(aligners):
        ax = axes[col_idx]
        for j, region in enumerate(GENOMICREGIONS):
            this_pcrt_avg = []
            for svtype in svtypes:
                this_pcrt_avg.append(np.mean(svtypes_pcrt[aligner][region][svtype]))

            if j == 0:
                ax.bar(xticks, this_pcrt_avg, color=REGIONCOLORS[region], width=bar_width, label=region)

            else:
                bottoms = []
                for k in range(0, j):
                    this_avg = []
                    for svtype in svtypes:
                        this_avg.append(np.mean(svtypes_pcrt[aligner][GENOMICREGIONS[k]][svtype]))
                    bottoms.append(this_avg)
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                ax.bar(xticks, this_pcrt_avg, bottom=bottom_sum, width=bar_width, color=REGIONCOLORS[region], label=region)

        ax.set_xticks(xticks)
        ax.set_xticklabels(svtypes, fontsize=13, rotation=60)

        ax.set_title(aligner, fontsize=13)

        if col_idx == 0:
            ax.set_ylabel('% of aligner unique SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(0, 80)
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 80, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig.tight_layout()
    plt.show()
