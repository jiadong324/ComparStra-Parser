#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/11/21
'''

import os
import pysam
import math
import pandas as pd

from Helpers.Annot import *
from Helpers.Constant import *
from Helpers.Functions import *


def compare_read_among_datasets():

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    for platform, datasets in DATASET_DICT.items():
        for caller in CALLERS:
            for aligner in ALIGNERS:
                ## Matches of all SVs between platforms
                merged_outdir = f'{WORKDIR}/read_dataset_repro/{platform}'
                if not os.path.exists(merged_outdir):
                    os.mkdir(merged_outdir)

                tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:

                    vcf_file = f'{WORKDIR}/{platform}/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligner} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                os.system(cmd)

                os.remove(tmp_file)

                get_dataset_compare_info(datasets, merged_out_vcf, merged_outdir, caller, aligner, platform, simple_reps, rmsk, sds)


def compare_assm_among_datasets():
    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')
    aligner = 'minimap2'

    plat_assembler = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['shasta', 'flye']}

    for platform, datasets in DATASET_DICT.items():

        for caller in ASMCALLERS:
            for assembler in plat_assembler[platform]:
                ## Matches of all SVs between platforms
                merged_outdir = f'{WORKDIR}/assm_dataset_repro/{platform}'
                if not os.path.exists(merged_outdir):
                    os.mkdir(merged_outdir)

                tmp_file = f'{merged_outdir}/{caller}.{assembler}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:
                    vcf_file = f'{WORKDIR}/{platform}/{aligner}_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{assembler}.jasmine.merged.vcf'
                print(f'Producing {caller} {assembler} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                os.system(cmd)

                os.remove(tmp_file)

                get_dataset_compare_info(datasets, merged_out_vcf, merged_outdir, caller, assembler, platform, simple_reps, rmsk, sds)

def get_dataset_compare_info(datasets, merged_vcf, compare_outdir, caller, tool, platform, simple_reps, rmsk, sds):

    matched_list = []
    unique_list = []

    extd_matched_list = []
    extd_unique_list = []

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

            if supp > 1:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_datasets = datasets[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_datasets, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_datasets, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

    # match_info_out = f'{workdir}/{sample}/{caller}/survivor/{sample}.{caller}.{aligner}.platform-concordant.{svtype}.tsv'
    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN', 'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}-concordant.info.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns =['ID_MATCH', 'TYPE_MATCH', 'DATASET', '#CHROM', 'POS', 'END', 'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}.unique.tsv', sep='\t', header=True, index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}-concordant.info.exTD.tsv', header=True, sep='\t',index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'DATASET', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}.unique.exTD.tsv', sep='\t', header=True, index=False)
