#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


import gzip
import pandas as pd

from Helpers.Functions import *
from Helpers.Annot import *
from Helpers.Constant import *

def process_assembly_calls():

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    exclude_dict = parse_exclude_regions(EXREGIONS)

    hifi_datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb']
    ont_datasets = ['ont_9kb', 'ont_19kb', 'ont_30kb']

    ont_assemblers = ['flye', 'shasta']
    hifi_assemblers = ['hifiasm', 'flye']

    pav_count_info = []
    pav_count_info_by_regions = []

    split_pav_vcf(hifi_datasets, hifi_assemblers, pav_count_info, pav_count_info_by_regions,
                  simple_reps, rmsk, sds, exclude_dict)

    split_pav_vcf(ont_datasets, ont_assemblers, pav_count_info, pav_count_info_by_regions,
                  simple_reps, rmsk, sds, exclude_dict)

    df_svcounts = pd.DataFrame(pav_count_info, columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num', 'cmrg_num', 'highconf_num'])
    df_svcounts.to_csv(f'{WORKDIR}/pav_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(pav_count_info_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{WORKDIR}/pav_sv_counts_region.tsv', sep='\t', header=True, index=False)

    svimasm_count_info = []
    svimasm_count_info_by_regions = []

    split_svimasm_vcf(hifi_datasets, hifi_assemblers, svimasm_count_info, svimasm_count_info_by_regions,
                      simple_reps, rmsk, sds, exclude_dict)

    split_svimasm_vcf(ont_datasets, ont_assemblers, svimasm_count_info,
                      svimasm_count_info_by_regions, simple_reps, rmsk, sds, exclude_dict)

    df_svcounts = pd.DataFrame(svimasm_count_info,columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num', 'cmrg_num', 'highconf_num'])
    df_svcounts.to_csv(f'{WORKDIR}/svimasm_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(svimasm_count_info_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{WORKDIR}/svimasm_sv_counts_region.tsv', sep='\t', header=True, index=False)


def split_pav_vcf(datasets, assemblers, sv_counts_info, sv_counts_by_regions,  simple_reps, rmsk, sds, exclude_dict):

    aligner = 'minimap2'
    svtypes = ['ins', 'del', 'inv']

    for assembler in assemblers:
        for dataset in datasets:
            current_dir = f'{WORKDIR}/HiFi'
            pav_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/pav_HG002.{assembler}.vcf.gz'
            if 'ont' in dataset:
                current_dir = f'{WORKDIR}/ONT'
                pav_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/pav_HG002.{assembler}.vcf.gz'

            total_num = 0

            pav_sv_vcf = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.vcf', 'w')
            bed_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.bed', 'w')

            svtype_by_region = {region: {SVTYPEMAP[svtype]: 0 for svtype in svtypes} for region in GENOMICREGIONS}

            ins_counter = 0
            del_counter = 0
            inv_counter = 0

            ins_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.ins.vcf', 'w')
            del_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.del.vcf', 'w')
            inv_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.pav.{assembler}.inv.vcf', 'w')

            with gzip.open(pav_vcf, 'rt') as f:
                for line in f:
                    if '#' in line:
                        print(line.strip(), file=ins_writer)
                        print(line.strip(), file=del_writer)
                        print(line.strip(), file=inv_writer)
                        print(line.strip(), file=pav_sv_vcf)
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

                    svtype = info_dict['SVTYPE']

                    if svtype == 'SNV':
                        continue

                    svlen = abs(int(info_dict['SVLEN']))

                    if svlen < 50 or svlen > MAXSIZE:
                        continue

                    this_chrom_exclude = exclude_dict[chrom]

                    end = start + svlen

                    if this_chrom_exclude.overlap(start, end):
                        continue

                    region, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simple_reps, rmsk, sds)

                    if svtype == 'INS':
                        ins_counter += 1

                        info_strs = ''
                        for key, val in info_dict.items():
                            info_strs += f'{key}={val};'
                        info_strs += f'END={end}'

                        str1 = '\t'.join(entries[0:7])
                        str2 = '\t'.join(entries[8:])
                        new_vcf_str = f'{str1}\t{info_strs}\t{str2}'

                        print(new_vcf_str, file=ins_writer)
                        print(new_vcf_str, file=pav_sv_vcf)

                        print(f'{chrom}\t{start}\t{end}\tINS\t{svlen}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=bed_writer)

                        svtype_by_region[region][svtype] += 1

                    elif svtype == 'DEL':
                        del_counter += 1
                        print(line.strip(), file=del_writer)
                        print(line.strip(), file=pav_sv_vcf)
                        print(f'{chrom}\t{start}\t{end}\tDEL\t{svlen}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=bed_writer)

                    elif svtype == 'INV':
                        inv_counter += 1
                        print(line.strip(), file=inv_writer)
                        print(line.strip(), file=pav_sv_vcf)
                        print(f'{chrom}\t{start}\t{end}\tINV\t{svlen}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=bed_writer)

                total_num += 1

            for region, sv_num in svtype_by_region.items():
                for svtype, count in sv_num.items():
                    sv_counts_by_regions.append((assembler, dataset, region, svtype, count))

            sv_counts_info.append(('pav', dataset, assembler, total_num, ins_counter, del_counter, inv_counter))

            print(f'Writing pav-minimap2-{dataset} #SVs: {total_num}')

            ins_writer.close()
            del_writer.close()
            inv_writer.close()
            pav_sv_vcf.close()

    df_svcounts = pd.DataFrame(sv_counts_info, columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num'])
    df_svcounts.to_csv(f'{WORKDIR}/pav_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(sv_counts_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{WORKDIR}/pav_sv_counts_region.tsv', sep='\t', header=True, index=False)

def split_svimasm_vcf(datasets, assemblers, sv_counts_info, sv_counts_by_regions, simple_reps, rmsk, sds, exclude_dict):

    aligner = 'minimap2'

    svtypes = ['ins', 'del', 'inv']

    for assembler in assemblers:
        for dataset in datasets:

            current_dir = f'{WORKDIR}/HiFi'
            input_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/variants.{assembler}.vcf'
            if 'ont' in dataset:
                current_dir = f'{WORKDIR}/ONT'
                input_vcf = f'{current_dir}/{aligner}_{dataset}/raw_calls/variants.{assembler}.vcf'

            print(f'Processing svimasm-{assembler}-{dataset} ...')

            total_num = 0
            ins_counter, del_counter, inv_counter = 0, 0, 0

            ins_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.ins.vcf', 'w')
            del_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.del.vcf', 'w')
            inv_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.inv.vcf', 'w')

            filtered_vcf_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.vcf', 'w')
            filtered_bed_writer = open(f'{current_dir}/{aligner}_{dataset}/filtered/HG002.svimasm.{assembler}.bed', 'w')

            svtype_by_region = {region: {SVTYPEMAP[svtype]: 0 for svtype in svtypes} for region in GENOMICREGIONS}

            for line in open(input_vcf, 'r'):
                if '#' in line:
                    print(line.strip(), file=del_writer)
                    print(line.strip(), file=ins_writer)
                    print(line.strip(), file=inv_writer)
                    print(line.strip(), file=filtered_vcf_writer)
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

                if svlen < 50 or svlen > MAXSIZE:
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

                elif svtype == 'DEL':
                    total_num += 1
                    del_counter += 1
                    print(line.strip(), file=del_writer)
                    print(line.strip(), file=filtered_vcf_writer)
                    svtype_by_region[region][svtype] += 1

                elif svtype == 'INV':
                    total_num += 1
                    inv_counter += 1
                    print(line.strip(), file=inv_writer)
                    print(line.strip(), file=filtered_vcf_writer)
                    svtype_by_region[region][svtype] += 1

                print(f'{chrom}\t{start}\t{end}\t{svtype}\t{svlen}\t{region}\t{rptype}\t{round(pcrt, 2)}', file=filtered_bed_writer)

            for region, sv_num in svtype_by_region.items():
                for svtype, count in sv_num.items():
                    sv_counts_by_regions.append((assembler, dataset, region, svtype, count))

            sv_counts_info.append(('svimasm', dataset, assembler, total_num, ins_counter, del_counter, inv_counter))

            del_writer.close()
            ins_writer.close()
            inv_writer.close()
            filtered_vcf_writer.close()
            filtered_bed_writer.close()


    df_svcounts = pd.DataFrame(sv_counts_info,columns=['caller', 'dataset', 'assembler', 'total', 'ins_num', 'del_num', 'inv_num'])
    df_svcounts.to_csv(f'{WORKDIR}/svimasm_sv_counts.tsv', sep='\t', header=True, index=False)

    df_svnum_regions = pd.DataFrame(sv_counts_by_regions, columns=['assembler', 'dataset', 'region', 'svtype', 'count'])
    df_svnum_regions.to_csv(f'{WORKDIR}/svimasm_sv_counts_region.tsv', sep='\t', header=True, index=False)

def merge_assm_insdel():

    caller_supp = 2

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    assembler_dict = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        assembler = assembler_dict[plat]
        for svtype in ['ins', 'del']:

            merged_caller_dir = f'{WORKDIR}/{plat}/assm_callers_merged'

            if not os.path.exists(merged_caller_dir):
                os.mkdir(merged_caller_dir)

            caller_vcf_path = f'{merged_caller_dir}/{dataset}_callers_sa2_vcf_path.txt'
            caller_vcf_path_writer = open(caller_vcf_path, 'w')

            for caller in ASMCALLERS:
                vcf_path = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{caller}.{assembler}.{svtype}.vcf'
                print(vcf_path, file=caller_vcf_path_writer)

            caller_vcf_path_writer.close()

            print(f' ==== Merge {svtype} on {dataset} ==== ')
            caller_sr2_merged_vcf = f'{merged_caller_dir}/{dataset}_{svtype}_callers_{assembler}_merged.vcf'

            cmd = f'{JASMINE} file_list={caller_vcf_path} out_file={caller_sr2_merged_vcf} max_dist=1000 --normalize_chrs --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1'

            os.system(cmd)
            os.remove(caller_vcf_path)
            obtain_confident_calls(caller_sr2_merged_vcf, merged_caller_dir, dataset, svtype, assembler, caller_supp, simple_reps, rmsk, sds)

def obtain_confident_calls(merged_vcf, workdir, dataset, svtype, assembler, supp_callers, simple_reps, rmsk, sds):
    suppvec_dict = {}
    merged_total = 0
    extd_supp_vec_dict = {}

    suppvec_info = {}
    extd_suppvec_info = {}

    extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged.extd.vcf', 'w')
    scs_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged.sc{supp_callers}.vcf', 'w')
    scs_extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged.sc{supp_callers}.extd.vcf', 'w')

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

            chrom, merged_id = entries[0], entries[2]
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

    supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged_suppinfo.txt', 'w')
    for supp, count in suppvec_dict.items():
        print(f'{supp}\t{count}', file=supp_writer)

    supp_writer.close()

    extd_supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged_suppinfo.extd.txt', 'w')
    for supp, count in extd_supp_vec_dict.items():
        print(f'{supp}\t{count}', file=extd_supp_writer)

    extd_supp_writer.close()