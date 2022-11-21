#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


from Reader.Assembly import *
from Reader.Read import *


def main():


    ####################################################
    ## Step1: Processing raw VCF files of each caller ##
    ####################################################

    process_read_calls()
    process_assembly_calls()

    ###########################################################
    ## Step2: Obtain high-confident insertions and deletions ##
    ###########################################################

    merge_reads_insdel()
    merge_assm_insdel()




if __name__ == '__main__':
    main()