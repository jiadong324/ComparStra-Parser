#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


import os
from intervaltree import IntervalTree
import scipy as sc

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%({v:d})'.format(p=pct,v=val)
    return my_autopct

def fmeasure(p, r):
    return 2*p*r / (p+r)


def intervene(bed_file_list, bed_file_names, output, bedtools=False):
    bed_file = ' '.join(bed_file_list)
    names = ','.join(bed_file_names)
    intervene_path = '/Users/apple/miniconda3/envs/dnatools/bin/intervene'
    cmd = f'{intervene_path} venn -i {bed_file} --names={names} -o {output}'
    if bedtools:
        cmd = f'{intervene_path} venn -i {bed_file} --names={names} -o {output} --bedtools-options f=0.5,r'
    os.system(cmd)


def parse_exclude_regions(exclude):
    exclude_dict = {}
    with open(exclude, 'r') as f:
        for line in f:
            entries = line.strip().split("\t")
            chrom = entries[0]

            start = int(entries[1])
            end = int(entries[2])
            if chrom not in exclude_dict:
                exclude_dict[chrom] = IntervalTree()
                exclude_dict[chrom][start:end] = (start, end)
            else:
                exclude_dict[chrom][start:end] = (start, end)
    return exclude_dict


def calculate_robust_score(c, n, m):
    return 1 - (2 * abs(c/n-c/m) / max(c/n, c/m))

def min_breakpoint_shift(begin_a, end_a, begin_b, end_b):
    if begin_b > end_a:
        bp_shift = begin_b - end_a

    else:
        if end_b < begin_a:
            bp_shift = begin_a - end_b
        else:
            bp_shift = min(abs(begin_a - begin_b), abs(end_a - end_b))

    return bp_shift

def max_breakpoint_shift(begin_a, end_a, begin_b, end_b):
    if begin_b > end_a:
        bp_shift = begin_b - end_a

    else:
        if end_b < begin_a:
            bp_shift = begin_a - end_b
        else:
            bp_shift = max(abs(begin_a - begin_b), abs(end_a - end_b))

    return bp_shift

def size_similarity(size_a, size_b):
    return round(min(size_a, size_b) / max(size_a, size_b), 3)

def reciprocal_overlap(begin_a, end_a, begin_b, end_b):
    overlap = min(end_a, end_b) - max(begin_a, begin_b)

    return round(min([overlap / (end_a - begin_a), overlap / (end_b - begin_b)]), 3)

def fmeasure_curve(f, p):
    return f * p / (2 * p - f)


def plot_fmeasures(ax, fstepsize=.1, stepsize=0.001):
    """Plots 10 fmeasure Curves into the current canvas."""
    p = sc.arange(0., 1., stepsize)[1:]
    for f in sc.arange(0., 1., fstepsize)[1:]:
        points = [(x, fmeasure_curve(f, x)) for x in p
                  if 0 < fmeasure_curve(f, x) <= 1.5]
        xs, ys = zip(*points)

        curve, = ax.plot(xs, ys, "--", color="gray", linewidth=1.7)  # , label=r"$f=%.1f$"%f) # exclude labels, for legend
        # bad hack:
        # gets the 10th last datapoint, from that goes a bit to the left, and a bit down
        ax.annotate(r"F1=%.1f" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.13, ys[-10] - 0.05), size=10)
        # ax.annotate(r"F1=%.1f" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.06, ys[-10] - 0.05), size=10)

def contains_gaps(chrom, start, end, ref):
    if start > end:
        start, end = end, start
    seq = ref.fetch(chrom, start - 2000, end + 2000)
    in_gap = False
    for base in seq:
        if base == 'N':
            in_gap = True
            break
    return in_gap

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
                return ('Tandem Repeats', subtype, overlap_pcrt)


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
                    return ('Tandem Repeats', subtype, overlap_pcrt)
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

def get_overlaps(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))