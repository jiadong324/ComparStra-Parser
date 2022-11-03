#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


import gzip

import matplotlib.pyplot as plt
import pandas as pd
import pysam
from statistics import mean
import seaborn as sns
import numpy as np

from Helpers.Functions import *
from Helpers.Annot import *
from Helpers.Constant import *

def process_assm_calls(workdir, max_size):

    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    exclude_dict = parse_exclude_regions(iMACEXCLUDE)
    cmrg_dict = parse_exclude_regions(CMRG)
    passed_dict = parse_exclude_regions(HIGHCONF)

    hifi_datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb']
    ont_datasets = ['ont_9kb', 'ont_19kb', 'ont_30kb']

    # hifi_datasets = ['hifi_18kb']
    # ont_datasets = ['ont_30kb']

    ont_assemblers = ['flye', 'shasta']
    hifi_assemblers = ['hifiasm', 'flye']

    pav_count_info = []
    pav_count_info_by_regions = []

    split_pav_vcf(workdir, hifi_datasets, hifi_assemblers, max_size, pav_count_info, pav_count_info_by_regions,
                  simple_reps, rmsk, sds, exclude_dict, cmrg_dict, passed_dict)

    split_pav_vcf(workdir, ont_datasets, ont_assemblers, max_size, pav_count_info, pav_count_info_by_regions,
                  simple_reps, rmsk, sds, exclude_dict, cmrg_dict, passed_dict)

    df_svcounts = pd.DataFrame(pav_count_info, columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num',
                                        'cmrg_num', 'highconf_num'])
    df_svcounts.to_csv(f'{workdir}/pav_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(pav_count_info_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{workdir}/pav_sv_counts_region.tsv', sep='\t', header=True, index=False)

    svimasm_count_info = []
    svimasm_count_info_by_regions = []

    split_svimasm_vcf(workdir, hifi_datasets, hifi_assemblers, max_size, svimasm_count_info, svimasm_count_info_by_regions,
                      simple_reps, rmsk, sds, exclude_dict, cmrg_dict, passed_dict)

    split_svimasm_vcf(workdir, ont_datasets, ont_assemblers, max_size, svimasm_count_info,
                      svimasm_count_info_by_regions, simple_reps, rmsk, sds, exclude_dict, cmrg_dict, passed_dict)

    df_svcounts = pd.DataFrame(svimasm_count_info,columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num', 'cmrg_num', 'highconf_num'])
    df_svcounts.to_csv(f'{workdir}/svimasm_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(svimasm_count_info_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{workdir}/svimasm_sv_counts_region.tsv', sep='\t', header=True, index=False)


def split_pav_vcf(workdir, datasets, assemblers, max_size, sv_counts_info, sv_counts_by_regions,  simple_reps, rmsk, sds, exclude_dict, cmrg_dict, passed_dict):

    aligner = 'minimap2'
    svtypes = ['ins', 'del', 'inv']

    for assembler in assemblers:
        for dataset in datasets:

            current_dir = f'{workdir}/HiFi'
            pav_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/pav_HG002.{assembler}.vcf.gz'
            if 'ont' in dataset:
                current_dir = f'{workdir}/ONT'
                pav_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/pav_HG002.{assembler}.vcf.gz'

            print(f'Processing pav-{assembler}-{dataset} ...')

            total_num = 0
            highconf_num = 0
            cmrg_num = 0

            pav_sv_vcf = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.vcf', 'w')
            bed_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.bed', 'w')

            # sv_in_cmrg =open(f'{current_dir}/{aligner}_{dataset}/filtered_in_cmrg/HG002.pav.{assembler}.vcf', 'w')
            # sv_in_passed = open(f'{current_dir}/{aligner}_{dataset}/filtered_in_highconf/HG002.pav.{assembler}.vcf', 'w')

            sv_header = 0

            svtype_num = {}
            svtype_by_region = {region: {SVTYPEMAP[svtype]: 0 for svtype in svtypes} for region in GENOMICREGIONS}
            for svtype in svtypes:
                counter = 0
                vcf_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.{svtype}.vcf', 'w')
                with gzip.open(pav_vcf, 'rt') as f:
                    for line in f:
                        if '#' in line:
                            print(line.strip(), file=vcf_writer)
                            if sv_header == 0:
                                print(line.strip(), file=pav_sv_vcf)
                                # print(line.strip(), file=sv_in_cmrg)
                                # print(line.strip(), file=sv_in_passed)
                            continue

                        entries = line.strip().split("\t")
                        chrom, start, sv_id = entries[0], int(entries[1]), entries[2]

                        if chrom not in AUTOSOMES:
                            continue

                        info_tokens = entries[7].split(";")
                        info_dict = {}

                        for token in info_tokens:
                            if "=" not in token:
                                continue
                            info_dict[token.split('=')[0]] = token.split('=')[1]

                        if info_dict['SVTYPE'] == 'SNV':
                            continue

                        svlen = abs(int(info_dict['SVLEN']))
                        if svlen < 50 or svlen > max_size:
                            continue

                        this_chrom_exclude = exclude_dict[chrom]

                        end = start + svlen

                        if this_chrom_exclude.overlap(start, end):
                            continue

                        if info_dict['SVTYPE'] == SVTYPEMAP[svtype]:
                            counter += 1

                            info_strs = ''
                            for key, val in info_dict.items():
                                info_strs += f'{key}={val};'
                            info_strs += f'END={end}'

                            str1 = '\t'.join(entries[0:7])
                            str2 = '\t'.join(entries[8:])
                            new_vcf_str = f'{str1}\t{info_strs}\t{str2}'

                            print(new_vcf_str, file=vcf_writer)
                            print(new_vcf_str, file=pav_sv_vcf)

                            region, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simple_reps, rmsk, sds)
                            print(f'{chrom}\t{start}\t{end}\t{SVTYPEMAP[svtype]}\t{svlen}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=bed_writer)

                            svtype_by_region[region][SVTYPEMAP[svtype]] += 1

                            # if cmrg_dict[chrom].overlap(start, end):
                            #     print(new_vcf_str, file=sv_in_cmrg)
                            #     cmrg_num += 1

                            # if passed_dict[chrom].overlap(start, end):
                            #     highconf_num += 1
                            #     print(new_vcf_str, file=sv_in_passed)


                vcf_writer.close()

                sv_header += 1

                total_num += counter

                svtype_num[svtype] = counter

            for region, sv_num in svtype_by_region.items():
                for svtype, count in sv_num.items():
                    sv_counts_by_regions.append((assembler, dataset, region, svtype, count))

            sv_counts_info.append(('pav', dataset, assembler, total_num, svtype_num['ins'], svtype_num['del'], svtype_num['inv'], cmrg_num, highconf_num))

            # bed_writer.close()
            pav_sv_vcf.close()
            # sv_in_cmrg.close()
            # sv_in_passed.close()

    df_svcounts = pd.DataFrame(sv_counts_info, columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num', 'cmrg_num', 'highconf_num'])
    df_svcounts.to_csv(f'{workdir}/pav_sv_counts.tsv', sep='\t', header=True, index=False)
    #
    df_svnum_regions = pd.DataFrame(sv_counts_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{workdir}/pav_sv_counts_region.tsv', sep='\t', header=True, index=False)

def split_svimasm_vcf(workdir, datasets, assemblers, max_size, sv_counts_info, sv_counts_by_regions, simple_reps, rmsk, sds, exclude_dict, cmrg_dict, passed_dict):

    aligner = 'minimap2'

    svtypes = ['ins', 'del', 'inv']

    for assembler in assemblers:
        for dataset in datasets:

            current_dir = f'{workdir}/HiFi'
            input_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/variants.{assembler}.vcf'
            if 'ont' in dataset:
                current_dir = f'{workdir}/ONT'
                input_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/variants.{assembler}.vcf'

            print(f'Processing svimasm-{assembler}-{dataset} ...')


            total_num = 0
            cmrg_num = 0
            highconf_num = 0
            ins_counter, del_counter, inv_counter = 0, 0, 0

            ins_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.ins.vcf', 'w')
            del_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.del.vcf', 'w')
            inv_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.inv.vcf', 'w')

            filtered_vcf_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.vcf', 'w')
            filtered_bed_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.bed', 'w')

            # sv_in_cmrg = open(f'{current_dir}/{aligner}_{dataset}/filtered_in_cmrg/5X/HG002.svimasm.{assembler}.vcf', 'w')
            # sv_in_passed = open(f'{current_dir}/{aligner}_{dataset}/filtered_in_highconf/5X/HG002.svimasm.{assembler}.vcf', 'w')

            svtype_by_region = {region: {SVTYPEMAP[svtype]: 0 for svtype in svtypes} for region in GENOMICREGIONS}

            for line in open(input_vcf, 'r'):
                if '#' in line:
                    print(line.strip(), file=ins_writer)
                    print(line.strip(), file=del_writer)
                    print(line.strip(), file=inv_writer)
                    print(line.strip(), file=filtered_vcf_writer)
                    # print(line.strip(), file=sv_in_passed)
                    # print(line.strip(), file=sv_in_cmrg)
                    continue

                entries = line.strip().split('\t')
                chrom, start, id = entries[0], int(entries[1]), entries[2]
                info_tokens = entries[7].split(';')
                info_dict = {}

                if chrom not in AUTOSOMES:
                    continue

                for token in info_tokens:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

                svtype = info_dict['SVTYPE']

                if svtype == 'BND':
                    continue

                end = int(info_dict['END'])
                svlen = abs(int(info_dict['SVLEN'])) if 'SVLEN' in info_dict else end - start

                if svlen < 50 or svlen > max_size:
                    continue

                end = start + svlen if svtype == 'INS' else int(info_dict['END'])

                this_chrom_exclude = exclude_dict[chrom]

                if this_chrom_exclude.overlap(start, end):
                    continue

                region, rptype, pcrt = annotate_sv_region(chrom, start, end, 50, simple_reps, rmsk, sds)

                if svtype == 'INS' or 'DUP' in svtype:
                    svtype = 'INS'
                    total_num += 1
                    ins_counter += 1
                    info_strs = ''
                    for key, val in info_dict.items():
                        if key == 'END':
                            continue
                        info_strs += f'{key}={val};'
                    info_strs += f'END={end}'

                    str1 = '\t'.join(entries[0:7])
                    str2 = '\t'.join(entries[8:])
                    new_vcf_str = f'{str1}\t{info_strs}\t{str2}'

                    print(new_vcf_str, file=ins_writer)
                    print(new_vcf_str, file=filtered_vcf_writer)
                    svtype_by_region[region][svtype] += 1

                    # if cmrg_dict[chrom].overlap(start, end):
                    #     cmrg_num += 1
                    #     print(new_vcf_str, file=sv_in_cmrg)

                    # if passed_dict[chrom].overlap(start, end):
                    #     highconf_num += 1
                    #     print(new_vcf_str, file=sv_in_passed)

                elif svtype == 'DEL':
                    total_num += 1
                    del_counter += 1
                    print(line.strip(), file=del_writer)
                    print(line.strip(), file=filtered_vcf_writer)
                    svtype_by_region[region][svtype] += 1

                    # if cmrg_dict[chrom].overlap(start, end):
                    #     cmrg_num += 1
                    #     print(line.strip(), file=sv_in_cmrg)

                    # if passed_dict[chrom].overlap(start, end):
                    #     highconf_num += 1
                    #     print(line.strip(), file=sv_in_passed)

                elif svtype == 'INV':
                    total_num += 1
                    inv_counter += 1
                    print(line.strip(), file=inv_writer)
                    print(line.strip(), file=filtered_vcf_writer)
                    svtype_by_region[region][svtype] += 1

                    # if cmrg_dict[chrom].overlap(start, end):
                    #     cmrg_num += 1
                    #     print(line.strip(), file=sv_in_cmrg)

                    # if passed_dict[chrom].overlap(start, end):
                    #     highconf_num += 1
                    #     print(line.strip(), file=sv_in_passed)

                print(f'{chrom}\t{start}\t{end}\t{svtype}\t{svlen}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=filtered_bed_writer)

            for region, sv_num in svtype_by_region.items():
                for svtype, count in sv_num.items():
                    sv_counts_by_regions.append((assembler, dataset, region, svtype, count))

            sv_counts_info.append(('svimasm', dataset, assembler, total_num, ins_counter, del_counter, inv_counter, cmrg_num, highconf_num))

            # bed_writer.close()
            ins_writer.close()
            del_writer.close()
            inv_writer.close()
            filtered_vcf_writer.close()
            filtered_bed_writer.close()

            # sv_in_cmrg.close()
            # sv_in_passed.close()


    df_svcounts = pd.DataFrame(sv_counts_info,columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num', 'cmrg_num', 'highconf_num'])
    df_svcounts.to_csv(f'{workdir}/svimasm_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(sv_counts_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{workdir}/svimasm_sv_counts_region.tsv', sep='\t', header=True, index=False)


def annotate_assm_highconf_svs(workdir, datasets, assemblers):

    simrep = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    aligner = 'minimap2'

    for assembler in assemblers:

        for dataset in datasets:
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in ASMCALLERS:
                print(f'Processing {aligner} {dataset} {caller} ....')
                vcf_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.{assembler}.vcf'
                annot_file = open(f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.{assembler}.annot.info.tsv', 'w')
                print('chrom\tstart\tend\tsvid\tsvlen\tsvtype\tregion\trptype\tpcrt', file=annot_file)

                with open(vcf_file, 'r') as f:
                    for line in f:
                        if '#' in line:
                            continue

                        entries = line.strip().split("\t")
                        chrom, start, sv_id = entries[0], int(entries[1]), entries[2]
                        info_tokens = entries[7].split(";")
                        info_dict = {}
                        for token in info_tokens:
                            if "=" not in token:
                                continue
                            info_dict[token.split("=")[0]] = token.split("=")[1].replace(">", "")
                        end = int(info_dict['END'])
                        sv_type = info_dict["SVTYPE"]
                        sv_len = end - start
                        if "SVLEN" in info_dict:
                            sv_len = abs(int(info_dict["SVLEN"]))
                        if "INS" in sv_type:
                            end = start + sv_len

                        region, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simrep, rmsk, sds)
                        print(f'{chrom}\t{start}\t{end}\t{sv_id}\t{sv_len}\t{sv_type}\t{region}\t{rptype}\t{pcrt}',
                              file=annot_file)

                annot_file.close()

def split_dipcall_vcf(workdir, datasets, aligners, max_size, simrep_path, rmsk_path, exclude_path, sd_path):
    simple_reps = pysam.Tabixfile(simrep_path, 'r')
    rmsk = pysam.Tabixfile(rmsk_path, 'r')
    sds = pysam.Tabixfile(sd_path, 'r')

    exclude_dict = parse_exclude_regions(exclude_path)

    sv_counts_info = []
    svlens = []
    for aligner in aligners:
        for dataset in datasets:

            print(f'Processing dipcall {aligner}-{dataset}...')

            input_vcf = f'{workdir}/{aligner}_{dataset}/HG002.dipcall.dip.vcf.gz'

            bed_dir = f'{workdir}/{aligner}_{dataset}/bed'

            if not os.path.exists(bed_dir):
                os.mkdir(bed_dir)

            bed_writer = open(f'{bed_dir}/{dataset}.dipcall.annot.bed', 'w')

            ins_counter, del_counter = 0, 0

            ins_writer = open(f'{workdir}/{aligner}_{dataset}/HG002.dipcall.ins.vcf', 'w')
            del_writer = open(f'{workdir}/{aligner}_{dataset}/HG002.dipcall.del.vcf', 'w')

            filtered_vcf_writer = open(f'{workdir}/{aligner}_{dataset}/HG002.dipcall.vcf', 'w')

            with gzip.open(input_vcf, 'rt') as f:
                for line in f:
                    if '#' in line:
                        print(line.strip(), file=ins_writer)
                        print(line.strip(), file=del_writer)
                        print(line.strip(), file=filtered_vcf_writer)
                        continue
                    entires = line.strip().split('\t')
                    chrom, start = entires[0], int(entires[1])
                    ref, alt = len(entires[3]), len(entires[4])
                    svlen = abs(alt - ref)
                    svlens.append(svlen)

                    if svlen < 50 or svlen > max_size:
                        continue

                    svtype = 'INS'
                    if alt - ref < 0:
                        svtype = 'DEL'

                    end = start + abs(alt - ref)
                    sv_id = f'{chrom}_{start}_{svtype}_{svlen}'

                    this_chrom_exclude = exclude_dict[chrom]

                    if this_chrom_exclude.overlap(start, end):
                        continue

                    infos = f'SVTYPE={svtype};SVLEN={svlen};END={end}'

                    outputs = f'{chrom}\t{start}\t{sv_id}\t{entires[3]}\t<{svtype}>\t{entires[5]}\t{infos}\t{entires[8]}\t{entires[9]}'
                    rptype, pcrt = rep_annotation(chrom, start, end, simple_reps, rmsk, sds)
                    print(f'{chrom}\t{start}\t{end}\t{svtype}\t{svlen}\t{rptype}\t{round(pcrt, 2)}', file=bed_writer)

                    if svtype == 'INS':
                        ins_counter += 1
                        print(outputs, file=ins_writer)
                        print(outputs, file=filtered_vcf_writer)
                    elif svtype == 'DEL':
                        del_counter += 1
                        print(outputs, file=del_writer)
                        print(outputs, file=filtered_vcf_writer)

            sv_counts_info.append(('dipcall', dataset, aligner, ins_counter + del_counter, ins_counter, del_counter, 0))
            ins_writer.close()
            del_writer.close()
            filtered_vcf_writer.close()
            bed_writer.close()

    df_svcounts = pd.DataFrame(sv_counts_info,columns=['caller', 'dataset', 'aligner', 'total', 'ins_num', 'del_num', 'inv_num'])
    df_svcounts.to_csv(f'{workdir}/dipcall_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svlens = pd.DataFrame(svlens, columns=['svlen'])

    fig, ax = plt.subplots(1,1, figsize=(6, 4))
    sns.histplot(df_svlens, x='svlen', log_scale=True, ax=ax)
    plt.show()