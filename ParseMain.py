#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''

from Reader.Assembly import *
from Reader.Read import *

from Com.Callers import compare_stra_callers
from Com.Datasets import compare_assm_among_datasets, compare_read_among_datasets
from Com.Aligners import compare_between_aligners
from Com.Assemblers import compare_between_assemblers
from Com.StraHQ import compare_stra_hq_insdel

def main():


    #############################################
    ## Processing raw VCF files of each caller ##
    ## Outputs are used to create Fig. 2a      ##
    #############################################

    process_read_calls()
    process_assembly_calls()

    ###############################################################
    ## Impact of dataset, aligner and assembler on each strategy ##
    ## Outputs are used to create Fig. 2b-2g                     ##
    ###############################################################

    compare_assm_among_datasets()
    compare_read_among_datasets()

    compare_between_aligners()
    compare_between_assemblers()


    ##############################################################################################
    ## Impact of aligner and assembler on SVs captured by read-based and assembly-based callers ##
    ## Outputs are used to create Fig. 3a-3c                                                    ##
    ##############################################################################################

    compare_stra_callers()


    ################################################################
    ## Obtain and compare high-confident insertions and deletions ##
    ## Outputs are used to create Fig. 3e-f and Fig. 4            ##
    ################################################################

    obtain_reads_hq_insdel()
    obtain_assm_hq_insdel()

    compare_stra_hq_insdel()



if __name__ == '__main__':
    main()