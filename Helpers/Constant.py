#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''

## Workpath
WORKDIR = '/home/jdlin/evaluate_caller/HG002/hs37d5'
EXREGIONS = '/home/jdlin/data/genome/grch37/grch37.exclude_regions_cen.bed'
HG19REF = '/home/jdlin/data/genome/grch37/hs37d5.fa'

SIMREP = '/home/jdlin/data/genome/grch37/simplerepeat.bed.gz'
RMSK = '/home/jdlin/data/genome/grch37/rmsk.bed.gz'
SD = '/home/jdlin/data/genome/grch37/seg_dup.bed.gz'

TRUEINSDEL = '/home/jdlin/data/genome/grch37/HG002_SVs_Tier1_v0.6.pass.bed'

AUTOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                  "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

SVTYPEMAP = {'ins': 'INS', 'del': 'DEL', 'others': 'Others', 'inv': 'INV', 'dup': 'DUP'}

GENOMICREGIONS = ['Tandem Repeats', 'Repeat Masked', 'Segment Dup', 'Simple Region']
CALLERS = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision']
ASMCALLERS = ['pav', 'svimasm']
ALIGNERS = ['minimap2', 'ngmlr', 'lra', 'winnowmap']
DATASETS = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']

TOOLMAP = {'pav': 'PAV', 'pbsv': 'pbsv', 'svim':'SVIM', 'cutesv': 'cuteSV', 'svision': 'SVision', 'sniffles':
    'Sniffles', 'nanovar': 'NanoVar', 'com': 'PAV', 'svimasm': 'SVIM-asm'}

MAXSIZE = 100000

## Tool path

JASMINE = ''
SAMTOOLS = ''
