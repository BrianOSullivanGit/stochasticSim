#!/bin/bash

# Spike in the required somatic variants to the phased pre-tumour BAMs then remap and combine both phased BAMs into the final output tumour BAM.
# Remap and combine phased normal BAMs into final output normal BAM.


# Usage:
# realignAndMerge.sh  <haplotype 1 BAM> \
#                         <haplotype 2 BAM>  \
#                           <reference> \
#                             <output BAM prefix> \
#
# example,
# sbatch -w node028 realignAndMerge.sh N_X1_HG00110.coding_exons_hg38_200x_76bp.bam N_X2_HG00110.coding_exons_hg38_200x_76bp.bam /workspace/user2/REF/GRCh38.d1.vd1.fa N_HG00110.coding_exons_hg38_200x_76bp

# Use standard tool versions linked in with this toolkit.
# User may use other versions in their path on their local machine if they wish.
# They will probably work, but as I havent tested it, I can not guarantee it.
# The path to the set of tools that have been tested is specified in a file
# tool.path that is located in ..../stochasticSim-x.y/bin.

if [ -z ${ALIGNER+xyz} ]; then ALIGNER=`which bwa`; fi
if [ -z ${SAMTOOLS+xyz} ]; then SAMTOOLS=`which samtools`; fi

# Setup the default number of cores we are using for this job (in samtools/bwa)
# Remember adding additional threads also means memory needs to be assigned also.
# More threads gets the job done faster but depending on your system you will run out of RAM if you don't get the balance right.
if [ -z ${SNGS_NUM_CORES+x} ]; then SNGS_NUM_CORES="4"; fi

hap1BamPrefix=`echo ${1} | sed -e 's/.*\///1' -e 's/\.bam$//1'`
hap2BamPrefix=`echo ${2} | sed -e 's/.*\///1' -e 's/\.bam$//1'`

# Reads coming from the simulation at this stage are 'perfectly' mapped.
# This is unrealistic.
# Remap BAMs from noth haplotypes and then merge them.

# Convert reads to fastq and remap with BWA/your aligner of choice..

# Convert to fastq.

# H1 set.
date | tr '\012' ':'
echo " Convert first haplotype BAM to fastq's and realign."
(${SAMTOOLS} collate -u -O ${1} || { echo -e "\n\033[7mSAMTOOLS COLLATE hap 1 failed! \033[0m"  1>&2;kill -HUP $$; }) | \
${SAMTOOLS} fastq -1 /tmp/temp${$}_1.fq -2 /tmp/temp${$}_2.fq -0 /dev/null -s /tmp/t${hap1BamPrefix}_S.fq -n || { echo -e "\n\033[7mSAMTOOLS FASTQ hap 1 failed! \033[0m"  1>&2;exit 1; }

date | tr '\012' ':'
echo " Realigning haplotype 1 using standard reference and convert to BAM."
(${ALIGNER} mem -t ${SNGS_NUM_CORES} -M ${3} /tmp/temp${$}_1.fq  /tmp/temp${$}_2.fq || { echo -e "\n\033[7mREALIGN hap 1 failed! \033[0m"  1>&2;kill -HUP $$; }) | \
(${SAMTOOLS} sort -@ ${SNGS_NUM_CORES} || { echo -e "\n\033[7mSAMTOOLS SORT hap 1 failed! \033[0m"  1>&2;exit 1; }) > ${hap1BamPrefix}.realn.bam

# Cleanup
rm ${hap1BamPrefix}.realn.sam /tmp/temp${$}_1.fq  /tmp/temp${$}_2.fq
# Index
${SAMTOOLS} index ${hap1BamPrefix}.realn.bam || { echo -e "\n\033[7mSAMTOOLS INDEX hap 1 failed! \033[0m"  1>&2;exit 1; }


# Now H2.
date | tr '\012' ':'
echo " Convert second haplotype BAM to fastq's and realign."
# mktemp temp${$}_1.fq temp${$}_2.fq

(${SAMTOOLS} collate -u -O ${2} || { echo -e "\n\033[7mSAMTOOLS COLLATE hap 2 failed! \033[0m"  1>&2;kill -HUP $$; }) | \
${SAMTOOLS} fastq -1 /tmp/temp${$}_1.fq -2 /tmp/temp${$}_2.fq -0 /dev/null -s ${hap2BamPrefix}_S.fq -n || { echo -e "\n\033[7mSAMTOOLS FASTQ hap 2 failed! \033[0m"  1>&2;exit 1; }


date | tr '\012' ':'
echo " Realigning haplotype 2 using standard reference and convert to BAM."
(${ALIGNER} mem -t ${SNGS_NUM_CORES} -M ${3} /tmp/temp${$}_1.fq  /tmp/temp${$}_2.fq || { echo -e "\n\033[7mREALIGN hap 2 failed! \033[0m"  1>&2;kill -HUP $$; }) | \
(${SAMTOOLS} sort -@ ${SNGS_NUM_CORES} || { echo -e "\n\033[7mSAMTOOLS SORT hap 2 failed! \033[0m"  1>&2;exit 1; })  > ${hap2BamPrefix}.realn.bam

# Cleanup
rm ${hap2BamPrefix}.realn.sam /tmp/temp${$}_1.fq  /tmp/temp${$}_2.fq
# Index
${SAMTOOLS} index ${hap2BamPrefix}.realn.bam || { echo -e "\n\033[7mSAMTOOLS INDEX hap 2 failed! \033[0m"  1>&2;exit 1; } 


date | tr '\012' ':'
echo " Merge both haplotypes into single BAM (${4}.realn.phased.bam)."
${SAMTOOLS} merge ${4}.realn.phased.bam ${hap1BamPrefix}.realn.bam ${hap2BamPrefix}.realn.bam || { echo -e "\n\033[7mSAMTOOLS MERGE failed! \033[0m"  1>&2;exit 1; } 

date | tr '\012' ':'
echo " Index the merged BAM."
${SAMTOOLS} index ${4}.realn.phased.bam || { echo -e "\n\033[7mSAMTOOLS INDEX failed! \033[0m"  1>&2;exit 1; }
