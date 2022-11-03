#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


import pandas as pd
import pysam
import gzip

from Helpers.Functions import *
from Helpers.Annot import *
from Helpers.Constant import *


def process_read_calls(datasets, aligners, max_size):
    exclude_dict = parse_exclude_regions(EXREGIONS)
    passed_regions = parse_exclude_regions(TRUEINSDEL)

    ref_file = pysam.FastaFile(HG19REF)

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    count_writer = open(f'{WORKDIR}/caller_sv_count.tsv', 'w')
    print('caller\tdataset\taligner\tall_num\tins_num\tdel_num\tinv_num\tdup_num\tbnd_num', file=count_writer)

    svtype_region = []
    svs_in_region = []

    for caller in CALLERS:
        for dataset in datasets:
            for aligner in aligners:
                current_dir = f'{WORKDIR}/{caller}/{aligner}_{dataset}'

                filtered_dir = f'{current_dir}/filtered'
                if not os.path.exists(filtered_dir):
                    os.mkdir(filtered_dir)

                vcf_path = f'{current_dir}/HG002.{caller}.s5.vcf'
                filtered_vcf_path = f'{filtered_dir}/HG002.{caller}.vcf'

                split_svtypes_prefix = f'HG002.{caller}'

                count_dict, svtype_by_region = split_svtypes_to_files(vcf_path, filtered_vcf_path, caller, -1, max_size, exclude_dict,
                                                    ref_file, filtered_dir, split_svtypes_prefix, simple_reps, rmsk, sds)

                all_svs, exbnds, ins_num, del_num, inv_num, dup_num, bnd_num = count_dict['Total'], count_dict['Exbnd'], \
                                                                         count_dict['INS'], count_dict['DEL'], \
                                                                         count_dict['INV'], count_dict['DUP'],count_dict['BND']

                print(f'{caller}\t{dataset}\t{aligner}\t{all_svs}\t{ins_num}\t{del_num}\t{inv_num}\t{dup_num}\t{bnd_num}', file=count_writer)

                print(f'Writing {caller}-{aligner}-{dataset} #SVs (exclude bnds): {all_svs}')

                for region, sv_num in svtype_by_region.items():
                    for svtype, count in sv_num.items():
                        svtype_region.append((caller, dataset, aligner, region, svtype, count))

                # cmrg_dir = f'{current_dir}/{aligner}_{dataset}/filtered_in_cmrg'
                # if not os.path.exists(cmrg_dir):
                #     os.mkdir(cmrg_dir)
                #
                # cmrg_svs, cmrg_svtypes = get_svs_in_cmrg(filtered_vcf_path, split_svtypes_prefix, cmrg_dir, cmrg_regions)
                # svs_in_region.append((caller, dataset, aligner, 'CMRG', cmrg_svs))


                high_conf_dir = f'{current_dir}/filtered_in_highconf'
                if not os.path.exists(high_conf_dir):
                    os.mkdir(high_conf_dir)

                passed_svs, passed_svtypes = get_svs_in_passed(filtered_vcf_path, split_svtypes_prefix, high_conf_dir, passed_regions)
                svs_in_region.append((caller, dataset, aligner, 'TRUE-INS/DEL', passed_svs))



    count_writer.close()

    df_svcounts = pd.DataFrame(svtype_region, columns=['caller', 'dataset', 'aligner', 'region', 'svtype', 'count'])
    df_svcounts.to_csv(f'{WORKDIR}/caller_sv_counts_region.tsv', sep='\t', header=True, index=False)

    df_svregions = pd.DataFrame(svs_in_region, columns=['caller', 'dataset', 'aligner', 'region_class', 'count'])
    df_svregions.to_csv(f'{WORKDIR}/caller_sv_in_ture_insdel.tsv', sep='\t', header=True, index=False)


def split_svtypes_to_files(input_vcf, noalt_vcf, caller, minsr, max_size, exclude_dict, ref_file, output_dir, output_prefix, simrep, rmsk, sds):

    re_tag_caller = ['sniffles', 'cutesv']
    filtered_sv_list = []
    sv_exbnd_list = []
    all_sv_list = []

    del_writer = open(f'{output_dir}/{output_prefix}.del.vcf', 'w')
    ins_writer = open(f'{output_dir}/{output_prefix}.ins.vcf', 'w')
    inv_writer = open(f'{output_dir}/{output_prefix}.inv.vcf', 'w')
    others_writer = open(f'{output_dir}/{output_prefix}.others.vcf', 'w')
    bnd_writer = open(f'{output_dir}/{output_prefix}.bnd.vcf', 'w')
    noalt_vcf_writer = open(noalt_vcf, 'w')

    exbnd_bed_writer = open(f'{output_dir}/{output_prefix}.exbnd.bed', 'w')

    # ins_counter = 0
    # del_counter = 0
    # inv_counter = 0
    # other_counter = 0
    bnd_counter = 0

    svtypes = ['INS', 'DEL', 'DUP', 'INV']
    svtype_count = {svtype: 0 for svtype in svtypes}
    svtype_by_region = {region: {svtype: 0 for svtype in svtypes} for region in GENOMICREGIONS}

    with open(input_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=del_writer)
                print(line.strip(), file=ins_writer)
                print(line.strip(), file=inv_writer)
                print(line.strip(), file=others_writer)
                print(line.strip(), file=bnd_writer)
                print(line.strip(), file=noalt_vcf_writer)
                continue

            entries = line.strip().split("\t")
            chrom, start, sv_id = entries[0], int(entries[1]), entries[2]

            if chrom not in AUTOSOMES:
                continue

            if entries[4][0] == '[' or entries[4][0] == ']':
                chr2 = entries[4][1: -2].split(':')[0]
                if chr2 not in AUTOSOMES:
                    continue
            if entries[4][-1] == ']' or entries[4][-1] == '[':
                chr2 = entries[4][2: -1].split(':')[0]
                if chr2 not in AUTOSOMES:
                    continue

            info_tokens = entries[7].split(";")
            info_dict = {}

            for token in info_tokens:
                if "=" not in token:
                    continue
                info_dict[token.split("=")[0]] = token.split("=")[1].replace(">", "")

            if minsr != -1:

                if caller in re_tag_caller and int(info_dict['RE']) < minsr:
                    continue

                if caller == 'svim' and int(info_dict['SUPPORT']) < minsr:
                    continue

            sv_type = info_dict["SVTYPE"]
            if 'DUP' in sv_type:
                sv_type = 'DUP'

            if sv_type in ['INS', 'DEL', 'INV', 'DUP']:
                end = int(info_dict["END"])
                sv_len = end - start
                if "SVLEN" in info_dict:
                    sv_len = abs(int(info_dict["SVLEN"]))

                if "INS" in sv_type:
                    end = start + sv_len

                if sv_len < 50 or sv_len >= max_size:
                    filtered_sv_list.append((chrom, start, end, sv_len))
                    continue

                if exclude_dict[chrom].overlap(start, end):
                    continue

                if contains_gaps(chrom, start, end, ref_file):
                    continue

                region, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simrep, rmsk, sds)

                svtype_by_region[region][sv_type] += 1
                svtype_count[sv_type] += 1

                if "INS" in sv_type:
                    info_strs = ''
                    for key, val in info_dict.items():
                        if key == 'END':
                            continue
                        info_strs += f'{key}={val};'
                    info_strs += f'END={end}'

                    str1 = '\t'.join(entries[0:7])
                    str2 = '\t'.join(entries[8:])
                    new_vcf_str = f'{str1}\t{info_strs}\t{str2}'

                    print(new_vcf_str, file=noalt_vcf_writer)


                else:
                    print(line.strip(), file=noalt_vcf_writer)

                print(f'{chrom}\t{start}\t{end}\t{sv_type}\t{sv_len}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=exbnd_bed_writer)

                if sv_type == 'INS' or sv_type == 'DUP':
                    print(line.strip(), file=ins_writer)

                elif sv_type == 'DEL':
                    print(line.strip(), file=del_writer)

                elif sv_type == 'INV':
                    print(line.strip(), file=inv_writer)

            elif sv_type == 'TRA' or sv_type == 'BND':
                bnd_counter += 1
                print(line.strip(), file=noalt_vcf_writer)
                print(line.strip(), file=bnd_writer)
                print(line.strip(), file=others_writer)
            else:
                print(line.strip(), file=others_writer)
                # print(f'{chrom}\t{start}\t{end}\t{sv_type}\t{sv_len}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=exbnd_bed_writer)


            all_sv_list.append((chrom, start, start + 1, sv_type, -1, sv_id))

    ins_writer.close()
    del_writer.close()
    inv_writer.close()
    bnd_writer.close()
    others_writer.close()
    noalt_vcf_writer.close()
    exbnd_bed_writer.close()

    count_dict = {'INS': svtype_count['INS'], 'DEL': svtype_count['DEL'], 'BND': bnd_counter, 'INV': svtype_count['INV'], 'DUP': svtype_count['DUP'],
                  'Total': len(all_sv_list), 'Exbnd': len(sv_exbnd_list)}

    return count_dict, svtype_by_region


def get_svs_in_passed(noalt_vcf, output_prefix, outdir, region_dict):

    svtypes = ['INS', 'DEL', 'DUP', 'INV']
    vcf_writer = open(f'{outdir}/{output_prefix}.vcf', 'w')
    counter = 0
    svtype_num = {svtype: 0 for svtype in svtypes}

    with open(noalt_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=vcf_writer)
                continue
            entries = line.strip().split("\t")
            chrom, start, sv_id = entries[0], int(entries[1]), entries[2]

            if chrom not in AUTOSOMES:
                continue

            if entries[4][0] == '[' or entries[4][0] == ']':
                chr2 = entries[4][1: -2].split(':')[0]
                if chr2 not in AUTOSOMES:
                    continue
            if entries[4][-1] == ']' or entries[4][-1] == '[':
                chr2 = entries[4][2: -1].split(':')[0]
                if chr2 not in AUTOSOMES:
                    continue

            info_tokens = entries[7].split(";")
            info_dict = {}

            for token in info_tokens:
                if "=" not in token:
                    continue
                info_dict[token.split("=")[0]] = token.split("=")[1].replace(">", "")

            sv_type = info_dict["SVTYPE"]
            if 'DUP' in sv_type:
                sv_type = 'DUP'

            if sv_type in svtypes:
                end = int(info_dict["END"])
                sv_len = end - start
                if "SVLEN" in info_dict:
                    sv_len = abs(int(info_dict["SVLEN"]))

                if "INS" in sv_type:
                    end = start + sv_len

                if region_dict[chrom].overlap(start, end):
                    counter += 1

                    if sv_type in svtype_num:
                        svtype_num[sv_type] += 1
                    else:
                        svtype_num[sv_type] = 1
                    print(line.strip(), file=vcf_writer)

    vcf_writer.close()
    return counter, svtype_num

def get_svs_in_cmrg(noalt_vcf, output_prefix, outdir, region_dict):


    svtypes = ['INS', 'DEL', 'DUP', 'INV']
    vcf_writer = open(f'{outdir}/{output_prefix}.vcf', 'w')
    counter = 0
    svtype_num = {svtype: 0 for svtype in svtypes}

    with open(noalt_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=vcf_writer)
                continue
            entries = line.strip().split("\t")
            chrom, start, sv_id = entries[0], int(entries[1]), entries[2]

            if chrom not in AUTOSOMES:
                continue

            if entries[4][0] == '[' or entries[4][0] == ']':
                chr2 = entries[4][1: -2].split(':')[0]
                if chr2 not in AUTOSOMES:
                    continue
            if entries[4][-1] == ']' or entries[4][-1] == '[':
                chr2 = entries[4][2: -1].split(':')[0]
                if chr2 not in AUTOSOMES:
                    continue

            info_tokens = entries[7].split(";")
            info_dict = {}

            for token in info_tokens:
                if "=" not in token:
                    continue
                info_dict[token.split("=")[0]] = token.split("=")[1].replace(">", "")

            sv_type = info_dict["SVTYPE"]
            if 'DUP' in sv_type:
                sv_type = 'DUP'

            if sv_type in svtypes:

                end = int(info_dict["END"])
                sv_len = end - start
                if "SVLEN" in info_dict:
                    sv_len = abs(int(info_dict["SVLEN"]))

                if "INS" in sv_type:
                    end = start + sv_len

                if region_dict[chrom].overlap(start, end):
                    counter += 1

                    if sv_type in svtype_num:
                        svtype_num[sv_type] += 1
                    else:
                        svtype_num[sv_type] = 1

                    print(line.strip(), file=vcf_writer)

    vcf_writer.close()

    return counter, svtype_num