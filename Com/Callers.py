#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/11/21
'''

import pickle

import pandas as pd
import pysam
import os
import math

from Helpers.Constant import *
from Helpers.Annot import *


def compare_stra_callers():
    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')
    regioned_svs_counts = []

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for assembler in plat_assemblers[plat]:
            for asm_method in ASMCALLERS:
                asm_calls = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{asm_method}.{assembler}.vcf'

                for caller in CALLERS:
                    for aligner in ALIGNERS:

                        print(f'Comparing {asm_method}-{assembler} to {caller}-{aligner} on {dataset} ...')

                        caller_calls = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                        compare_outdir = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/comstra_callers'

                        if not os.path.exists(compare_outdir):
                            os.mkdir(compare_outdir)

                        tmp_file = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.txt'
                        tmp_file_writer = open(tmp_file, 'w')
                        print(asm_calls, file=tmp_file_writer)
                        print(caller_calls, file=tmp_file_writer)

                        tmp_file_writer.close()

                        merged_vcf = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{dataset}.jasmine.merged.vcf'
                        cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                        os.system(cmd)

                        os.remove(tmp_file)
                        supp_info = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{dataset}.jasmine.suppinfo.tsv'
                        write_suppvec_info(merged_vcf, supp_info)

                        svs_by_regions = get_caller_compare_info([asm_method, caller], merged_vcf, compare_outdir, caller, asm_method, aligner, assembler, readset, simple_reps, rmsk, sds)

                        for region_label, counts in svs_by_regions.items():
                            regioned_svs_counts.append((caller, asm_method, dataset, aligner, assembler, region_label, counts[0], counts[1], counts[2]))

    df_regioned_svs = pd.DataFrame(regioned_svs_counts, columns=['caller', 'asm_method', 'dataset', 'aligner', 'assembler', 'region', 'assm_unique','intersects', 'align_unique'])
    df_regioned_svs.to_csv(f'{WORKDIR}/strategy_compare_byregions.tsv', header=True, sep='\t', index=False)


def get_caller_compare_info(callers, merged_vcf, compare_outdir, read_caller, asm_caller, aligner, assembler, dataset, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    extd_matched_list = []
    extd_unique_list = []

    svs_by_regions = {repeat: [0, 0, 0] for repeat in GENOMICREGIONS}


    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue

            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if supp == len(callers):
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_caller = callers[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_caller, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))
                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_caller, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))

                if supp_vec == '01':
                    svs_by_regions[region_label][2] += 1
                if supp_vec == '10':
                    svs_by_regions[region_label][0] += 1


    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.tsv', sep='\t', header=True, index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.exTD.tsv',header=True, sep='\t', index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.exTD.tsv', sep='\t', header=True, index=False)

    return svs_by_regions