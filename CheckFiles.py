#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/24

'''
import logging
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Constant import *

def check_files():


    print('\n===== Check VCF file completeness =====')

    if not os.path.exists(HG19REF):
        logging.error('Cannot locate reference file, please provide the correct path!')

    if not os.path.exists(SIMREP):
        logging.error('Cannot locate simplerepeat.bed.gz file, please provide the correct path!')

    if not os.path.exists(RMSK):
        logging.error('Cannot locate rmsk.bed.gz file, please provide the correct path!')

    if not os.path.exists(SD):
        logging.error('Cannot locate seg_dup.bed.gz file, please provide the correct path!')

    if not os.path.exists(EXREGIONS):
        logging.error('Cannot locate grch37.exclude_regions_cen.bed file, please provide the correct path!')

    if not os.path.exists(f'{WORKDIR}/CMRGs'):
        logging.error(f'Cannot find {WORKDIR}/CMRGs')

    if not os.path.exists(f'{WORKDIR}/CMRGs/plat_compare'):
        logging.error(f'Cannot find {WORKDIR}/CMRGs/plat_compare')

    if not os.path.exists(f'{WORKDIR}/CMRGs/stra_compare'):
        logging.error(f'Cannot find {WORKDIR}/CMRGs/stra_compare')

    if not os.path.exists(f'{WORKDIR}/truvari'):
        logging.error(f'Cannot find {WORKDIR}/truvari')

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    counter = 0
    for dataset in DATASETS:
        platform = 'HiFi'
        if 'ont' in dataset:
            platform = 'ONT'

        if not os.path.exists(f'{WORKDIR}/{platform}/ExTD_hq_trio_validation'):
            logging.error(f'Cannot find {WORKDIR}/{platform}/ExTD_hq_trio_validation')

        for caller in READCALLERS:
            for aligner in READALIGNERS:
                current_dir = f'{WORKDIR}/{platform}/{aligner}_{dataset}'
                vcf_path = f'{current_dir}/raw_calls/HG002.{caller}.s5.vcf'
                counter += 1
                if not os.path.exists(vcf_path):
                    logging.error(f'Missed VCF file: {aligner} {dataset} {TOOLMAP[caller]}')

        for aligner in ['minimap2', 'lra']:
            for assembler in plat_assemblers[platform]:
                pav_vcf = f'{WORKDIR}/{platform}/{aligner}_{dataset}/raw_calls/pav_HG002.{assembler}.vcf.gz'
                svimasm_vcf = f'{WORKDIR}/{platform}/{aligner}_{dataset}/raw_calls/variants.{assembler}.vcf'
                counter += 2
                if aligner == 'lra' and dataset in ['hifi_10kb', 'hifi_15kb']:
                    continue

                if not os.path.exists(pav_vcf):
                    logging.error(f'Missed VCF file: {aligner} {dataset} PAV')


                if not os.path.exists(svimasm_vcf):
                    logging.error(f'Missed VCF file: {aligner} {dataset} SVIM-asm')



    print(f'All {counter} VCF files and other extra files are ready for analysis!!')


if __name__ == '__main__':
    check_files()
