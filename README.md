
## Amis

Evaluating the impact of assemblers, aligners, sequencer and read length on read-based and assembly-based SV detection.

## Analysis workflow

### Overview 

The major parts involved in the comparison are listed below:
1. Using six long read datasets ([details](https://github.com/jiadong324/ComparStra-Parser/wiki)), we assessed the impact of dataset, aligner and assembler on the detection variability.
2. On each dataset, 20 read-based callsets and four assembly-based callsets were compared to assess the impact of aligner and assembler.
3. Based on the analysis of 2, we build high-confident insertions and deletions (insdel) callsets of read and assembly. The high-confident insdel callsets are then compared.
4. Benchmarking 20 read-based and eight assembly-based detection piplines with well curate SVs of HG002 released by GIAB.

### Project files at Zenodo

Please check the [wiki](https://github.com/jiadong324/ComparStra-Parser/wiki) page for more details about working directory structure, [SV detection](https://github.com/jiadong324/ComparStra-Parser/wiki/SV-detection) and [benchmarking](https://github.com/jiadong324/ComparStra-Parser/wiki/HG002-benchmarking).

There are 12 separate zip files to download from Zenodo ([![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7856049.svg)](https://doi.org/10.5281/zenodo.7856049)
).


1. ```lra_ont_9kb.zip```
2. ```lra_ont_19kb.zip```
3. ```lra_ont_30kb.zip```
4. ```minimap2_ont_9kb.zip```
5. ```minimap2_ont_19kb.zip```
6. ```minimap2_ont_30kb.zip```
7. ```other_files_under_ONT.zip```
8. ```HiFi dataset calls.zip```
9. ```CMRGs.zip```
10. ```truvari.zip```
11. ```hg19_ref.zip```
12. ```hg19_repeats.zip```

Unzip ```hg19_ref.zip``` to get Hg19 reference genome hs37d5.fa.

Unzip ```hg19_repeats.zip``` to get files listed below that are used in the analysis. Please refer [CAMPHOR](https://github.com/afujimoto/CAMPHOR) to get more details for processing the repeat files.

1) Simple repeat file including STR and VNTR (simplerepeat.bed.gz). 
2) Segmental duplication file (seg_dup.bed.gz).
3) Repeat masker file including LINE, SINE and etc (rmsk.bed.gz).
4) Hg19 excluded regions (grch37.exclude_regions_cen.bed).

The above files from 1-10 are also available at OneDrive ([reproduce_data.zip](https://stuxjtueducn-my.sharepoint.com/:f:/g/personal/pengjia_stu_xjtu_edu_cn/EnORz16WyTpEtWGWCT8IziUB8Rtm_x-iUGm2UQ4u6bV5Sw?e=MhMGxB)).

Unzip ```reproduce_data.zip``` and you will get all files under the working directory ```reproduce_data```.

### Analysis environment

#### Required tools and packages

```
## Tools
Jasmine=1.1.4
Samtools=1.9

## Python packages
python=3.6
pandas=1.1.5
numpy=1.19.5
seaborn=0.11.1
pysam=0.15.3
matplotlib_venn=0.11.7
intervaltree=3.1.0
```

#### Create environment for data analysis

```
## Create a python environment
conda create -n py36 python=3.6
conda activate py36

## Install required packages
pip install seaborn==0.11.1
pip install matplotlib-venn==0.11.9
pip install pysam==0.15.3
pip install intervaltree==3.1.0

## Install Jasmine
conda config --add channels bioconda
conda config --add channels conda-forge
conda install jasminesv
```


### Reproducing results

#### Setup variables

Please assign the absolute path to the following variables in ```./Helpers/Constant.py```

```
WORKDIR = '/path/to/reproduce_data'
FIGDIR = '/path/to/reproduce_data/Figures'

HG19REF = '/path/to/hs37d5.fa'
EXREGIONS = '/path/to/hg19_repeats/grch37.exclude_regions_cen.bed'
SIMREP = '/path/to/hg19_repeats/simplerepeat.bed.gz'
RMSK = '/path/to/hg19_repeats/rmsk.bed.gz'
SD = '/path/to/hg19_repeats/seg_dup.bed.gz'

SAMTOOLS = '/path/to/samtools'
JASMINE = '/path/to/jasmine'
```

**NOTE:** Please run the scripts by the order listed below. 

#### Figure 2
```
## Figure 2a
python ./Figure2/Figure2a.py
## Figure 2b and 2c
python ./Figure2/Figure2bc.py
## Figure 2d, 2e, 2f and 2g
python ./Figure2/Figure2defg.py
```

#### Figure 3
```
## Figure 3a, 3b and 3c
python ./Figure3/Figure3abc.py
## Figure 3d, 3e and 3f
python ./Figure3/Figure3def.py
```

#### Figure 4
```
## Figure 4a, 4b, 4c and 4d
python ./Figure4/Figure4.py
```

#### Figure 5
```
## Figure 5a, 5b, 5c, 5d, 5e and 5f
python ./Figure5/Figure5.py
```

#### Extended Data Figures

```
## Extended Data Fig 1
python ./SuppFig/FigS1.py

## Extended Data Fig 2
python ./SuppFig/FigS2.py

## Extended Data Fig 3
python ./SuppFig/FigS3.py

## Extended Data Fig 4
python ./SuppFig/FigS4.py

## Extended Data Fig 5
python ./SuppFig/FigS5.py

## Extended Data Fig 6
python ./SuppFig/FigS6.py

## Extended Data Fig 7
python ./SuppFig/FigS7.py
```

