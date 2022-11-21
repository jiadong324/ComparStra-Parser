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


def process_read_calls():

    exclude_dict = parse_exclude_regions(EXREGIONS)
    ref_file = pysam.FastaFile(HG19REF)

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    count_writer = open(f'{WORKDIR}/caller_sv_count.tsv', 'w')
    print('caller\tdataset\taligner\tall_num\tins_num\tdel_num\tinv_num\tdup_num\tbnd_num', file=count_writer)

    svtype_region = []
    svs_in_region = []

    for caller in CALLERS:
        for dataset in DATASETS:
            for aligner in ALIGNERS:
                current_dir = f'{WORKDIR}/{caller}/{aligner}_{dataset}'

                filtered_dir = f'{current_dir}/filtered'
                if not os.path.exists(filtered_dir):
                    os.mkdir(filtered_dir)

                vcf_path = f'{current_dir}/HG002.{caller}.s5.vcf'
                filtered_vcf_path = f'{filtered_dir}/HG002.{caller}.vcf'

                split_svtypes_prefix = f'HG002.{caller}'

                count_dict, svtype_by_region = split_svtypes_to_files(vcf_path, filtered_vcf_path, caller, -1, exclude_dict,
                                                    ref_file, filtered_dir, split_svtypes_prefix, simple_reps, rmsk, sds)

                all_svs, exbnds, ins_num, del_num, inv_num, dup_num, bnd_num = count_dict['Total'], count_dict['Exbnd'], \
                                                                         count_dict['INS'], count_dict['DEL'], \
                                                                         count_dict['INV'], count_dict['DUP'],count_dict['BND']

                print(f'{caller}\t{dataset}\t{aligner}\t{all_svs}\t{ins_num}\t{del_num}\t{inv_num}\t{dup_num}\t{bnd_num}', file=count_writer)

                print(f'Writing {caller}-{aligner}-{dataset} #SVs (exclude bnds): {all_svs}')

                for region, sv_num in svtype_by_region.items():
                    for svtype, count in sv_num.items():
                        svtype_region.append((caller, dataset, aligner, region, svtype, count))

    count_writer.close()

    df_svcounts = pd.DataFrame(svtype_region, columns=['caller', 'dataset', 'aligner', 'region', 'svtype', 'count'])
    df_svcounts.to_csv(f'{WORKDIR}/caller_sv_counts_region.tsv', sep='\t', header=True, index=False)

    df_svregions = pd.DataFrame(svs_in_region, columns=['caller', 'dataset', 'aligner', 'region_class', 'count'])
    df_svregions.to_csv(f'{WORKDIR}/caller_sv_in_ture_insdel.tsv', sep='\t', header=True, index=False)


def split_svtypes_to_files(input_vcf, noalt_vcf, caller, minsr, exclude_dict, ref_file, output_dir, output_prefix, simrep, rmsk, sds):

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

                if sv_len < 50 or sv_len >= MAXSIZE:
                    filtered_sv_list.append((chrom, start, end, sv_len))
                    continue

                if exclude_dict[chrom].overlap(start, end):
                    continue

                if contains_gaps(chrom, start, end, ref_file):
                    continue

                region, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simrep, rmsk, sds)

                svtype_by_region[region][sv_type] += 1
                svtype_count[sv_type] += 1

                if "INS" in sv_type or "DUP" in sv_type:
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
                    print(new_vcf_str, file=ins_writer)


                    print(f'{chrom}\t{start}\t{end}\t{sv_type}\t{sv_len}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=exbnd_bed_writer)

                elif sv_type == 'DEL':
                    print(line.strip(), file=del_writer)
                    print(line.strip(), file=noalt_vcf_writer)

                elif sv_type == 'INV':
                    print(line.strip(), file=inv_writer)
                    print(line.strip(), file=noalt_vcf_writer)

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


def merge_reads_insdel():
    caller_supp = 5

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for svtype in ['ins', 'del']:
            sr4_caller_vcfs = f'{WORKDIR}/{plat}/{dataset}_callers_sr4_vcf_path.txt'
            sr4_caller_writer = open(sr4_caller_vcfs, 'w')

            for caller in CALLERS:
                vcf_path = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{caller}.{svtype}.vcf'

                print(vcf_path, file=sr4_caller_writer)

            sr4_caller_writer.close()

            merged_caller_dir = f'{WORKDIR}/{plat}/read_callers_merged'

            print(f' ==== Merge {svtype} on {dataset} ==== ')
            caller_sr4_merged_vcf = f'{merged_caller_dir}/{dataset}_{svtype}_callers_minimap2_merged.vcf'
            cmd = f'{JASMINE} file_list={sr4_caller_vcfs} out_file={caller_sr4_merged_vcf} max_dist=1000 --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1'

            os.system(cmd)
            os.remove(sr4_caller_vcfs)
            obtain_confident_calls(caller_sr4_merged_vcf, merged_caller_dir, dataset, svtype, caller_supp, simple_reps,
                                   rmsk, sds)


def obtain_confident_calls(merged_vcf, workdir, dataset, svtype, supp_callers, simple_reps, rmsk, sds):
    suppvec_dict = {}
    merged_total = 0
    extd_supp_vec_dict = {}

    suppvec_info = {}
    extd_suppvec_info = {}

    extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged.extd.vcf', 'w')
    scs_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged.sc{supp_callers}.vcf', 'w')
    scs_extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged.sc{supp_callers}.extd.vcf', 'w')

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=extd_merged_vcf)
                print(line.strip(), file=scs_merged_vcf)
                print(line.strip(), file=scs_extd_merged_vcf)
                continue

            entries = line.strip().split('\t')
            info_tokens = entries[7].split(";")
            info_dict = {}
            merged_total += 1

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            supp_vec = info_dict['SUPP_VEC']
            callers = []
            for i, val in enumerate(supp_vec):
                if val == '1':
                    callers.append(TOOLMAP[CALLERS[i]])

            merged_id = entries[2]

            suppvec_info[merged_id] = ','.join(callers)

            if supp_vec in suppvec_dict:
                suppvec_dict[supp_vec] += 1
            else:
                suppvec_dict[supp_vec] = 1

            merged_type = info_dict['SVTYPE']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if region_label != 'Tandem Repeats':
                print(line.strip(), file=extd_merged_vcf)

                extd_suppvec_info[merged_id] = ','.join(callers)

                if supp_vec in extd_supp_vec_dict:
                    extd_supp_vec_dict[supp_vec] += 1
                else:
                    extd_supp_vec_dict[supp_vec] = 1

            if supp >= supp_callers:
                print(line.strip(), file=scs_merged_vcf)
                if region_label != 'Tandem Repeats':
                    print(line.strip(), file=scs_extd_merged_vcf)

    supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged_suppinfo.txt', 'w')
    for supp, count in suppvec_dict.items():
        print(f'{supp}\t{count}', file=supp_writer)

    supp_writer.close()

    extd_supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged_suppinfo.extd.txt', 'w')
    for supp, count in extd_supp_vec_dict.items():
        print(f'{supp}\t{count}', file=extd_supp_writer)

    extd_supp_writer.close()