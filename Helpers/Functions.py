#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


import os
import pandas as pd
import vcf
from intervaltree import IntervalTree
from math import pi
import matplotlib.pylab as plt
import scipy as sc

from Helpers.Constant import *

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

def get_survivor_supp(merged_vcf):

    supp_dict = {}
    merged_total = 0
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            try:
                info_tokens = entries[7].split(";")
                info_dict = {}
                merged_total += 1

                for token in info_tokens:
                    if '=' in token:
                        info_dict[token.split('=')[0]] = token.split('=')[1]

                supp = int(info_dict['SUPP'])
                if supp in supp_dict:
                    supp_dict[supp] += 1
                else:
                    supp_dict[supp] = 1
            except Exception:
                print(line)

    return supp_dict, merged_total

def get_survivor_suppvec(merged_vcf):

    suppvec_dict = {}
    merged_total = 0
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_tokens = entries[7].split(";")
            info_dict = {}
            merged_total += 1

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = info_dict['SUPP_VEC']
            if supp in suppvec_dict:
                suppvec_dict[supp] += 1
            else:
                suppvec_dict[supp] = 1

    return suppvec_dict, merged_total

def read_suppvec_info(supp_info):
    suppvec_dict = {}
    total = 0
    with open(supp_info, 'r') as f:
        for line in f:
            entries = line.strip().split('\t')
            if entries[0] == 'Total':
                total = int(entries[1])
            else:
                suppvec_dict[entries[0]] = int(entries[1])
    return suppvec_dict, total