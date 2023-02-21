#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/11/21
'''

import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import vcf
import os, math

import upsetplot
from matplotlib.lines import Line2D
import pysam
from statistics import mean
from matplotlib.patches import Patch


from Helpers.Constant import *
from Helpers.Annot import annotate_sv_region



def compare_between_aligners():

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in CALLERS:
            outdir = f'{WORKDIR}/{plat}/{dataset}_aligner_repro'
            if not os.path.exists(outdir):
                os.mkdir(outdir)

            tmp_file = f'{outdir}/{caller}.{dataset}.txt'
            tmp_file_writer = open(tmp_file, 'w')

            for aligner in ALIGNERS:
                vcf_file = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                print(f'{vcf_file}', file=tmp_file_writer)
            tmp_file_writer.close()

            merged_out_vcf = f'{outdir}/{caller}.{dataset}.jasmine.merged.vcf'
            print(f'Producing {caller} {dataset} merged calls ...')
            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1'
            os.system(cmd)

            print(f'Annotating {caller} {dataset} ...')
            get_aligner_compare_info(ALIGNERS, merged_out_vcf, outdir, caller, simple_reps, rmsk, sds)

            os.remove(tmp_file)


def get_aligner_compare_info(aligners, merged_vcf, compare_outdir, caller, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    extd_matched_list = []
    extd_unique_list = []

    svs_by_regions = {rep: [0 for i in range(len(aligners))] for rep in GENOMICREGIONS}

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

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)

            if supp > 1:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_aligner = aligners[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_aligner, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_aligner, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

                if supp_vec == '01':
                    svs_by_regions[region_label][2] += 1
                if supp_vec == '10':
                    svs_by_regions[region_label][0] += 1


    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{caller}.aligner-concordant.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{caller}.aligner-unique.tsv', sep='\t', header=True, index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{caller}.aligner-concordant.exTD.tsv', header=True, sep='\t', index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{caller}.aligner-unique.exTD.tsv', sep='\t', header=True, index=False)