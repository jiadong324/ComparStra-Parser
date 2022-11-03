#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/9/26

'''


from Reader.Assm import *
from Reader.Read import *


from Helpers.Constant import *


def main():

    max_size = 100000
    datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']
    aligners = ['minimap2', 'ngmlr', 'lra', 'winnowmap']
    hifi_datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb']
    ont_datasets = ['ont_9kb', 'ont_19kb', 'ont_30kb']

    ont_assemblers = ['shasta', 'flye']
    hifi_assemblers = ['hifiasm', 'flye']

    '''
        Process raw calls of each caller
    '''

    process_read_calls(datasets, aligners, max_size)
    # process_assm_calls(WORKDIR, max_size)

    # assm_binned_tree, assm_binned_list = create_bins()

    # annotate_read_highconf_svs(iMACDIR, datasets, aligners)
    # annotate_assm_highconf_svs(iMACDIR, ont_datasets, ont_assemblers)
    # annotate_assm_highconf_svs(iMACDIR, hifi_datasets, hifi_assemblers)


    '''
        Merge FNs/FPs at CMRGs
        1. Merge all caller FNs of each dataset (stra.dataset.fn.merged.vcf).
        2. Obtain all-caller-fn of each dataset (xxx.HiFi.all-caller-fn.merged.vcf, xxx.ONT.all-caller-fn.merged.vcf).
        3. Merge all-caller-fn among all datasets (stra.HiFi.fn.merged.vcf, stra.ONT.fn.merged.vcf)
    '''
    # merge_cmrg(f'{iMAC}/HG002/CMRGs/truvari_results', dataset_dict, 'fp')

if __name__ == '__main__':
    main()