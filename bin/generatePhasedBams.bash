#!/bin/bash

# Generate phased reads from personalised diploid target (ie., genome which only includes target region,
# plus any off-target flanking regions if specified)...
# A set of phased normal and pre-tumour BAMs are created as output.

# generatePhasedBams.sh <haplotype 1 reference> <haplotype 2 reference> <read length> <-f|-c> <fold coverage | number of reads> <mean fragment length> <user suffix>

# Example use:
# generatePhasedBams.sh X1_HG00110.coding_exons_hg38.fa X2_HG00110.coding_exons_hg38.fa 76 143000000 180 "200x_76bp"
# generatePhasedBams.sh X1_HG00110.coding_exons_hg38.fa X2_HG00110.coding_exons_hg38.fa 76 250250000 180 "350x_76bp"

# Outputs normal and pre-tumour BAMs

# Use standard tool versions linked in with this toolkit.
# User may use other versions in their path on their local machine if they wish.
# They will probably work, but as I havent tested it, I can not guarantee it.
# The path to the set of tools that have been tested is specified in a file
# tool.path that is located in ..../stochasticSim-x.y/bin.
if [ -z ${SAMTOOLS+xyz} ]; then SAMTOOLS=`which samtools`; fi
if [ -z ${READ_SIMULATOR+xyz} ]; then READ_SIMULATOR=`which art_illumina"`; fi

hap1RefPrefix=`echo ${1} | sed -e 's/.*\///1' -e 's/\.fa$//1'`
hap2RefPrefix=`echo ${2} | sed -e 's/.*\///1' -e 's/\.fa$//1'`

# Simulate reads, sam format.
# Remember if you change READ_SIMULATOR above you obviously need to
# update command line switches below accordingly.
date | tr '\012' ':'
echo " Simulating phased reads for pre-tumour sample.."
${READ_SIMULATOR} -na -sam -i ${1} -p -l ${3} ${4} -m ${5} -s 10 -o "T_${hap1RefPrefix}_${6}"
#${READ_SIMULATOR} -na -sam -i ${1} -p -l ${3} -c ${4} -m ${5} -s 10 -o "T_${hap1RefPrefix}_${6}"

# Remove the paired end fastq files, we won't need those.
rm T_${hap1RefPrefix}_${6}[12].fq

${READ_SIMULATOR} -na -sam -i ${2} -p -l ${3} ${4} -m ${5} -s 10 -o "T_${hap2RefPrefix}_${6}"
#${READ_SIMULATOR} -na -sam -i ${2} -p -l ${3} -c ${4} -m ${5} -s 10 -o "T_${hap2RefPrefix}_${6}"

# Again, remove redundant paired end fastq files.
rm T_${hap2RefPrefix}_${6}[12].fq

date | tr '\012' ':'
echo " Prefix QNAME in tumour reads with haplotype identifier (haplotype 1 / T_X1_HG00110)."
cat T_${hap1RefPrefix}_${6}.sam | sed -e "s/^[^@]/T_X1_&/1" > T_${hap1RefPrefix}_${6}.tmp.sam

date | tr '\012' ':'
echo " Prefix QNAME in tumour reads with haplotype identifier (haplotype 2 / T_X2_HG00110)."
cat T_${hap2RefPrefix}_${6}.sam | sed -e "s/^[^@]/T_X2_&/1" > T_${hap2RefPrefix}_${6}.tmp.sam


# Sort, index and convert to (perfectly mapped) phased pre-tumour BAMs.
date | tr '\012' ':'
echo " Sort, index and convert to T_${hap1RefPrefix}_${6}.bam."
${SAMTOOLS} sort T_${hap1RefPrefix}_${6}.tmp.sam > T_${hap1RefPrefix}_${6}.bam
${SAMTOOLS} index T_${hap1RefPrefix}_${6}.bam

date | tr '\012' ':'
echo " Sort, index and convert to T_${hap2RefPrefix}_${6}.bam."
${SAMTOOLS} sort T_${hap2RefPrefix}_${6}.tmp.sam > T_${hap2RefPrefix}_${6}.bam
${SAMTOOLS} index T_${hap2RefPrefix}_${6}.bam

# Cleanup
rm [TN]_*_${6}*.sam

##### NORMAL

date | tr '\012' ':'
echo " Simulating phased reads from normal sample.."
${READ_SIMULATOR} -na -sam -i ${1} -p -l ${3} ${4} -m ${5} -s 10 -o "N_${hap1RefPrefix}_${6}"
#${READ_SIMULATOR} -na -sam -i ${1} -p -l ${3} -c ${4} -m ${5} -s 10 -o "N_${hap1RefPrefix}_${6}"

# Remove the paired end fastq files, we won't need those.
rm N_${hap1RefPrefix}_${6}[12].fq

${READ_SIMULATOR} -na -sam -i ${2} -p -l ${3} ${4} -m ${5} -s 10 -o "N_${hap2RefPrefix}_${6}"
#${READ_SIMULATOR} -na -sam -i ${2} -p -l ${3} -c ${4} -m ${5} -s 10 -o "N_${hap2RefPrefix}_${6}"

# Remove the paired end fastq files, we won't need those.
rm N_${hap2RefPrefix}_${6}[12].fq


date | tr '\012' ':'
echo " Prefix QNAME in normal reads with haplotype identifier (haplotype 1 / N_X1_HG00110)."
cat N_${hap1RefPrefix}_${6}.sam | sed -e "s/^[^@]/N_X1_&/1" > N_${hap1RefPrefix}_${6}.tmp.sam

date | tr '\012' ':'
echo " Prefix QNAME in normal reads with haplotype identifier (haplotype 2 / N_X2_HG00110)."
cat N_${hap2RefPrefix}_${6}.sam | sed -e "s/^[^@]/N_X2_&/1" > N_${hap2RefPrefix}_${6}.tmp.sam


# Sort, index and convert to (perfectly mapped) phased normal BAMs.
date | tr '\012' ':'
echo " Sort, index and convert to N_${hap1RefPrefix}_${6}.bam."
${SAMTOOLS} sort N_${hap1RefPrefix}_${6}.tmp.sam > N_${hap1RefPrefix}_${6}.bam
${SAMTOOLS} index N_${hap1RefPrefix}_${6}.bam

date | tr '\012' ':'
echo " Sort, index and convert to N_${hap2RefPrefix}_${6}.bam."
${SAMTOOLS} sort N_${hap2RefPrefix}_${6}.tmp.sam > N_${hap2RefPrefix}_${6}.bam
${SAMTOOLS} index N_${hap2RefPrefix}_${6}.bam


# Cleanup
rm [TN]_*_${6}*.sam
