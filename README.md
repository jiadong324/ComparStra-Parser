
## Amis

Evaluating the impact of assemblers, aligners, sequencer and read length on read-based and assembly-based SV detection.

## Datasets

| Dataset | Description | Source | 
| --- | ----------- | ---- |
| HiFi-10kb | 10kb average read length | GIAB [FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/)|
| HiFi-15kb | 15kb average read length | GIAB [FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/)|
| HiFi-18kb | 18kb average read length| SRA Acession IDs: SRR18239004, SRR18239005, SRR18239006 and SRR18239007|
| ONT-9kb | 9kb average read length | GIAB [FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/) |
| ONT-19kb | 19kb average read length | SRA Acession ID: SRR18363746|
| ONT-30kb | 30kb average read length | SRA Acession ID: SRR18363752 and SRR18363753|

## Analysis workflow


### Overview 

The major parts involved in the comparison are listed below:
1. For each strategy, we assessed the impact of dataset, aligner and assembler on the detection variability.
2. On each dataset, 20 read-based callsets and four assembly-based callsets were compared to assess the impact of aligner and assembler.
3. Based on the analysis of 2, we build high-confident insertions and deletions (insdel) callsets of read and assembly. The high-confident insdel callsets are then compared.
4. Benchmarking 20 read-based and eight assembly-based detection piplines with well curate SVs of HG002 released by GIAB.

### Reads alignment

#### minimap2
```
## HiFi
minimap2 -a -H -k 19 -O 5,56 -E 4,1 -A 2 -B 5 -z 400,50 -r 2000 -g 5000 -Y --MD 

## ONT
minimap2 -a -z 600,200 -x map-ont -Y --MD 
```

#### ngmlr

```
## HiFi
ngmlr -x pacbio

## ONT
ngmlr -x ont
```

#### winnowmap

```
## HiFi
winnowmap -W ./repetitive_k15.txt -ax map-pb -Y --MD

## ONT
winnowmap -W ./repetitive_k15.txt -ax map-ont -Y --MD
```

#### lra

```
## HiFi
zcat $fastq | lra align -CCS $ref /dev/stdin -t $thread -p s --printMD -SkipH 

## ONT
zcat $fastq | lra align -ONT $ref /dev/stdin -t $thread -p s --printMD -SkipH

```

### Sequence assembly and alignment

#### flye

```
## HiFi
flye --pacbio-hifi

## ONT
flye --nano-raw
```

#### shasta

```
./shasta-Linux-0.8.0 --config Nanopore-Oct2021
```

#### hifiasm

```
hifiasm -o output.tag -l1 input.fqs input.fqs_untag
```

#### Assembly alignment

The alignment parameters are used in PAV.
```
minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 --secondary=no -O 5,56 -E 4,1 -B 5 -a --eqx -Y
```

### SV detection


For all read-based callers, the minimum number of support read is set to 1, 2, 5 and 5 for 5X, 10X, 20X and 35X coverage data, respectively

#### pbsv

Minimum number of support read is set with parameters ```-A``` and ```-O```.
```
pbsv discover
pbsv call -m 50 -A 5 -O 5 -S 0
```

#### SVIM

Minimum number of support read is set with bcftools command ``` bcftools view -i "SUPPORT >= X" ```.
```
svim alignment --cluster_max_distance 1.4 --min_sv_size 50 

## Minimum number of support read 2 for 10X and 20X coverage data
bcftools view -i "SUPPORT >= 5" variants.vcf > HG002.svim.s5.vcf
```

#### cuteSV

Minimum number of support read is set with parameter ```-s```.
```
## HiFi
cuteSV --min_size 50 -t 6 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -s 5

## ONT
cuteSV --min_size 50 -t 6 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3  -s 5 
```

#### SVision

Minimum number of support read is set with parameter ```-s```.
```
SVision -n HG002 -s 5
```

#### Sniffles

Minimum number of support read is set with parameter ```--minsupport```.
```
sniffles --minsvlen 50 --minsupport 5
```

#### PAV

```
snakemake -s Snakefile  -j 28  -k --ri >sublog 2>&1 &
```

#### SVIM-asm

```
svim-asm diploid --tandem_duplications_as_insertions --interspersed_duplications_as_insertions 
```

### HG002 benchmarking

#### Post-processing

Example of create compressed and indexed VCF file for further evaluation.
```
cat HG002.pav.flye.vcf| awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | bgzip -c > ./pav.flye.sorted.vcf.gz
tabix ./pav.flye.sorted.vcf.gz
```
**NOTE:** 
1. For PAV output, the header for SVLEN should be changed to ```<Number=1,Type=Integer>```
1. For SVision output, the ```Covered``` in the filter column should be replaced with ```PASS```

#### Truvari evaluation

**NOTE:** We do not consider genotype accuracy for benchmarking and the Truvari version was v3.0.0.

For SVs at true INS/DEL regions ([download link](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/))

```
truvari bench -f ./hs37d5.fa -b ./HG002_SVs_Tier1_v0.6.vcf.gz --includebed ./HG002_SVs_Tier1_v0.6.bed --passonly --giabreport -r 1000 -p 0.00 -c /pav.flye.sorted.vcf.gz -o pav_flye
```

For SV at CMRGs ([download link](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant/))

```
truvari bench -f ./hs37d5.fa -b ./HG002_GRCh37_CMRG_SV_v1.00.vcf.gz --includebed ./HG002_GRCh37_CMRG_SV_v1.00.bed --passonly --giabreport -r 1000 -p 0.00 -c /pav.flye.sorted.vcf.gz -o pav_flye
```