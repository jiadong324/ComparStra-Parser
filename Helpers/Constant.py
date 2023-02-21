#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''

## PKL-workpath
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
ASSMBLERCOLOR = {'flye': '#005f86', 'hifiasm': '#c06600', 'shasta': '#737a00'}

STRACOLORS = {'Assembly': '#97a7bc', 'Read': '#f4a578'}
SHIFTCOLORDICT = {'0': '#acaed6', '0,10':'#e9b9c0', '10,50': '#dba038', '>10':'#949e80'}

TOOLMARKERS = {'pav': 'P', 'pbsv': 'X', 'svim': 'd', 'cutesv': 'o', 'sniffles': 's', 'svision': 'H',
               'svimasm': 'p'}
TOOLMARKERS2 = {'PAV': 'P', 'pbsv': 'X', 'SVIM': 'd', 'cuteSV': 'o', 'Sniffles': 's', 'SVision': 'H',
               'SVIM-asm': 'p'}

TOOLCOLORS = {'svimasm': '#c078aa', 'pbsv': '#d95f02', 'svim': '#7570b3', 'cutesv': '#d51dad', 'svision': '#008066',
              'sniffles': '#386cb0', 'pav': '#c88f1c', 'dipcall': '#e41a1c'}
TOOLCOLORS2 = {'SVIM-asm': '#c078aa', 'pbsv': '#d95f02', 'SVIM': '#7570b3', 'cuteSV': '#d51dad', 'SVision': '#008066',
              'Sniffles': '#386cb0', 'PAV': '#c88f1c'}

TOOLMAP = {'pav': 'PAV', 'pbsv': 'pbsv', 'svim':'SVIM', 'cutesv': 'cuteSV', 'svision': 'SVision', 'sniffles':
    'Sniffles', 'nanovar': 'NanoVar', 'com': 'PAV', 'svimasm': 'SVIM-asm', 'dipcall': 'dipcall'}

PLATMAP = {'hifi_10kb': 'HiFi-10kb', 'hifi_15kb': 'HiFi-15kb',
            'hifi_18kb': 'HiFi-18kb', 'ont_9kb': 'ONT-9kb', 'ont_19kb': 'ONT-19kb', 'ont_30kb': 'ONT-30kb',}

PLATCOLORS = {'HiFi': '#F2984E', 'ONT': '#8EAB8F'}
ALIGNERCOLOR = {'minimap2': '#ce1181', 'ngmlr': '#1b9e77', 'lra': '#d1664d', 'winnowmap': '#85abd3'}
REGIONCOLORS = {'Tandem Repeats': '#7599c1', 'Repeat Masked': '#65936a', 'Segment Dup': '#916e8c', 'Simple Region': '#d1cbae'}

SVTYPECOLORS = {'DEL': '#557b8a', 'DUP': '#00a330', 'INS': '#e06c41', 'INV': '#f2ed4f', 'BND': '#c265e7', 'Others': '#CDB699', 'DUP:TANDEM': '#CDB699'}
