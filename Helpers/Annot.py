#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''

import pysam
import vcf

from Helpers.Constant import *


def annotate_sv_region(chrom, start, end, pcrt_thresh, simreps_tabix, rmsk_tabix, sd_tabix):
    if start > end:
        start, end = end, start
    size = end - start + 1
    annotations = []

    if 'chr' not in chrom:
        chrom = f'chr{chrom}'

    for simrep in simreps_tabix.fetch(chrom, start, end):
        entries = simrep.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)

        if overlap_size > 0:
            motif = rp_info.split(',')[-1]
            overlap_pcrt = min(overlap_size / size * 100, 100)
            subtype = 'VNTR' if len(motif) >= 7 else 'STR'
            if overlap_pcrt >= pcrt_thresh:
                annotations.append(('Tandem Repeats', subtype, overlap_pcrt))


    for rmsk in rmsk_tabix.fetch(chrom, start, end):
        entries = rmsk.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[4]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            rptype = rp_info.split(',')[11]
            if overlap_pcrt >= pcrt_thresh:
                if rptype == 'Simple_repeat':
                    motif = rptype[1: -2]
                    subtype = 'VNTR' if len(motif) >= 7 else 'STR'
                    annotations.append(('Tandem Repeats', subtype, overlap_pcrt))
                    continue
                annotations.append(('Repeat Masked', rptype, overlap_pcrt))

    for sd in sd_tabix.fetch(chrom, start, end):
        entries = sd.strip().split('\t')
        sd_start, sd_end, sd_mate_coord = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, sd_start, sd_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            annotations.append(('Segment Dup', 'SegDup', overlap_pcrt))

    if len(annotations) == 0:
        return ('Simple Region', 'None', 0)

    sorted_annotations = sorted(annotations, key=lambda x:x[1], reverse=True)

    return sorted_annotations[0]


def rep_annotation(chrom, start, end, simreps_tabix, rmsk_tabix, sd_tabix):
    if start > end:
        start, end = end, start
    size = end - start + 1
    annotations = []

    if 'chr' not in chrom:
        chrom = f'chr{chrom}'

    for simrep in simreps_tabix.fetch(chrom, start, end):
        entries = simrep.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)

        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            motif = rp_info.split(',')[-1]
            if len(motif) >= 7:
                annotations.append(('VNTR', overlap_pcrt))
            else:
                annotations.append(('STR', overlap_pcrt))

    for rmsk in rmsk_tabix.fetch(chrom, start, end):
        entries = rmsk.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[4]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            rptype = rp_info.split(',')[11]
            if rptype == 'Simple_repeat':
                continue
            annotations.append((rptype, overlap_pcrt))


    for sd in sd_tabix.fetch(chrom, start, end):
        entries = sd.strip().split('\t')
        sd_start, sd_end, sd_mate_coord = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, sd_start, sd_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            annotations.append(('SegDup', overlap_pcrt))


    if len(annotations) == 0:
        return ('None', 0)

    sorted_annotations = sorted(annotations, key=lambda x:x[1], reverse=True)

    return sorted_annotations[0]

def get_overlaps(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))



def rep_annotate_uniques(outdir, asm_caller, caller, svtype, merged_vcf_path, simple_reps, rmsk, sds):
    merged_vcf = vcf.Reader(open(merged_vcf_path, 'r'))

    asm_uniques = open(f'{outdir}/{asm_caller}_uniques.{svtype}.annot.tsv', 'w')
    align_uniques = open(f'{outdir}/{caller}_uniques.{svtype}.annot.tsv', 'w')

    for rec in merged_vcf:
        supp_vec = rec.INFO['SUPP_VEC']
        chrom, start = rec.CHROM, int(rec.POS)
        end, svlen = int(rec.INFO['END']), int(rec.INFO['SVLEN'])

        if start == end:
            end += svlen

        if supp_vec == '10':
            pav_id = rec.samples[0].data[7]
            rptype, pcrt = rep_annotation(chrom, start, end, simple_reps, rmsk, sds)
            print(f'{chrom}\t{start}\t{pav_id}\t{svlen}\t{svtype}\t{rptype}\t{round(pcrt, 2)}\tminimap2', file=asm_uniques)

        elif supp_vec == '01':

            rptype, pcrt = rep_annotation(chrom, start, end, simple_reps, rmsk, sds)
            print(f'{chrom}\t{start}\t{rec.ID}\t{svlen}\t{svtype}\t{rptype}\t{round(pcrt, 2)}\tminimap2',file=align_uniques)

    asm_uniques.close()
    align_uniques.close()