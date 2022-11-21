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

from Helpers.Constant import *
from Helpers.Functions import *

def compare_stra_hq_insdel():

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')
    assembler_dict = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    for plat, dataset in DATASET_DICT.items():
        assembler = assembler_dict[plat]

        for svtype in ['ins', 'del']:
            stra_compar_dir = f'{WORKDIR}/{plat}/stra_compare_hq/WGS'
            if not os.path.exists(stra_compar_dir):
                os.mkdir(stra_compar_dir)

            stra_compare_vcf_path = f'{stra_compar_dir}/stra_compare_path.txt'
            stra_compare_vcf_path_writer = open(stra_compare_vcf_path, 'w')
            print(f'{WORKDIR}/{plat}/assm_callers_merged/{dataset}_{svtype}_callers_{assembler}_merged.sc2.vcf', file=stra_compare_vcf_path_writer)
            print(f'{WORKDIR}/{plat}/read_callers_merged/{dataset}_{svtype}_callers_minimap2_merged.sc5.vcf', file=stra_compare_vcf_path_writer)

            stra_compare_vcf_path_writer.close()

            stra_compare_merged_vcf = f'{stra_compar_dir}/{dataset}_{svtype}_assm_read_merged.WGS.vcf'

            cmd = f'{JASMINE} file_list={stra_compare_vcf_path} out_file={stra_compare_merged_vcf} max_dist=1000 --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1'
            os.system(cmd)
            os.remove(stra_compare_vcf_path)

            get_hq_compare_info(dataset, svtype, stra_compare_merged_vcf, stra_compar_dir, 'WGS', simple_reps, rmsk, sds)

            extd_stra_compar_dir = f'{WORKDIR}/{plat}/stra_compare_hq/ExTD'
            if not os.path.exists(extd_stra_compar_dir):
                os.mkdir(extd_stra_compar_dir)

            extd_stra_compare_vcf_path = f'{extd_stra_compar_dir}/extd_stra_compare_path.txt'
            extd_stra_compare_vcf_path_writer = open(extd_stra_compare_vcf_path, 'w')
            print(f'{WORKDIR}/{plat}/assm_callers_merged/{dataset}_{svtype}_callers_{assembler}_merged.sc2.extd.vcf',
                  file=extd_stra_compare_vcf_path_writer)
            print(f'{WORKDIR}/{plat}/read_callers_merged/{dataset}_{svtype}_callers_minimap2_merged.sc5.extd.vcf',
                  file=extd_stra_compare_vcf_path_writer)

            extd_stra_compare_vcf_path_writer.close()
            extd_stra_compare_merged_vcf = f'{extd_stra_compar_dir}/{dataset}_{svtype}_assm_read_merged.ExTD.vcf'
            cmd = f'{JASMINE} file_list={extd_stra_compare_vcf_path} out_file={extd_stra_compare_merged_vcf} max_dist=1000 --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1'

            os.system(cmd)
            os.remove(extd_stra_compare_vcf_path)

            get_hq_compare_info(dataset, svtype, extd_stra_compare_merged_vcf, extd_stra_compar_dir, 'ExTD', simple_reps, rmsk, sds)


def get_hq_compare_info(dataset, svtype, merged_vcf, compare_outdir, label, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    svs_by_regions = {repeat: [0, 0, 0] for repeat in GENOMICREGIONS}
    supp_vec_dict = {}

    if not os.path.exists(compare_outdir):
        os.mkdir(compare_outdir)

    read_unique_bed = open(f'{compare_outdir}/{dataset}.{svtype}.read-uniques.{label}.bed', 'w')
    assm_unique_bed = open(f'{compare_outdir}/{dataset}.{svtype}.assm-uniques.{label}.bed', 'w')
    assm_unique_vcf = open(f'{compare_outdir}/{dataset}.{svtype}.assm-uniques.{label}.vcf', 'w')
    assm_unique_counter = 0
    read_unique_counter = 0

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=assm_unique_vcf)
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
            sv_length = info_dict['SVLEN']

            if supp_vec in supp_vec_dict:
                supp_vec_dict[supp_vec] += 1
            else:
                supp_vec_dict[supp_vec] = 1


            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if supp == 2:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_list.append((merged_id, merged_type, supp_vec, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))
                if supp_vec == '01':
                    read_unique_counter += 1
                    svs_by_regions[region_label][2] += 1

                    if svtype == 'ins':
                        print(f'{entries[0]}\t{entries[1]}\t{entries[1]}\t{merged_id}\tINS_{sv_length}', file=read_unique_bed)
                    else:
                        print(f'{entries[0]}\t{entries[1]}\t{end}\t{merged_id}\tDEL', file=read_unique_bed)

                if supp_vec == '10':
                    assm_unique_counter += 1
                    svs_by_regions[region_label][0] += 1
                    print(line.strip(), file=assm_unique_vcf)
                    if svtype == 'ins':
                        print(f'{entries[0]}\t{entries[1]}\t{entries[1]}\t{merged_id}\tINS_{sv_length}\t', file=assm_unique_bed)
                    else:
                        print(f'{entries[0]}\t{entries[1]}\t{end}\t{merged_id}\tDEL', file=assm_unique_bed)

    assm_unique_bed.close()
    read_unique_bed.close()
    assm_unique_vcf.close()

    print(f'Read specific {svtype}: {read_unique_counter}')
    print(f'Assembly specific {svtype}: {assm_unique_counter}')

    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{dataset}.{svtype}.concordants.info.{label}.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{dataset}.{svtype}.uniques.info.{label}.tsv', sep='\t', header=True, index=False)

    suppvec_writer = open(f'{compare_outdir}/{dataset}.{svtype}.suppvec_info.{label}.tsv', 'w')
    for supp_vec, count in supp_vec_dict.items():
        print(f'{supp_vec}\t{count}', file=suppvec_writer)
    suppvec_writer.close()
