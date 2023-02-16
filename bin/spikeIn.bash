#!/bin/bash

# Spike in the required somatic variants to the phased pre-tumour BAMs.
# After this is complete you may remap and combine phased tumour BAMs into final output tumour BAM.

# spikeIn.bash  <haplotype 1 pre-tumour BAM> \
#                         <haplotype 2 pre-tumour BAM>  \
#                           <haplotype 1 reference> \
#                             <haplotype 2 reference> \
#                               <haplotype 1 tumour somatics list> \
#                                 <haplotype 2 tumour somatics list>

# (for example)
# spikeIn.sh   T_X1_HG00110.coding_exons_hg38_30x_76bp.bam T_X2_HG00110.coding_exons_hg38_30x_76bp.bam X1_HG00110.coding_exons_hg38.fa X2_HG00110.coding_exons_hg38.fa X1_HG00110_coding_exons_hg38.merged.spike X2_HG00110_coding_exons_hg38.merged.spike


# see pointmassSpikeInCfg.sh for how config (the .spike file) is created.


# TODO!!!! somaticSim below takes a seed if you want it... Its hardcoded to 434, change this to an option if someone wants it.

# Use standard tool versions linked in with this toolkit.
# User may use other versions in their path on their local machine if they wish.
# They will probably work, but as I havent tested it, I can not guarantee it.
# The path to the set of tools that have been tested is specified in a file
# tool.path that is located in ..../stochasticSim-x.y/bin.

if [ -z ${SAMTOOLS+xyz} ]; then SAMTOOLS=`which samtools`; fi
if [ -z ${STOCHASTIC_SPIKE+xyz} ]; then STOCHASTIC_SPIKE=`which stochasticSpike`; fi

hap1PreTumourBamPrefix=`echo ${1} | sed -e 's/.*\///1' -e 's/\.bam$//1'`
#echo $hap1PreTumourBamPrefix

hap2PreTumourBamPrefix=`echo ${2} | sed -e 's/.*\///1' -e 's/\.bam$//1'`
#echo $hap2PreTumourBamPrefix


# Spike in these, each set into its corresponding haplotype set of reads.
date | tr '\012' ':'
echo " Spike into haplotype 1 ${1}"
${STOCHASTIC_SPIKE} ${1} ${3} ${5} 434 ${hap1PreTumourBamPrefix}.spike.sam || { echo -e "\n\033[7mSPIKE IN hap 1 failed! \033[0m";exit 1; }
mv truth.vcf ${hap1PreTumourBamPrefix}.truth.vcf

date | tr '\012' ':'
echo " Spike into haplotype 2 ${2}"
${STOCHASTIC_SPIKE} ${2} ${4} ${6} 434 ${hap2PreTumourBamPrefix}.spike.sam || { echo -e "\n\033[7mSPIKE IN hap 2 failed! \033[0m";exit 1; }
mv truth.vcf ${hap2PreTumourBamPrefix}.truth.vcf

date | tr '\012' ':'
echo " Sort, index and convert to ${hap1PreTumourBamPrefix}.spike.bam."
${SAMTOOLS} sort ${hap1PreTumourBamPrefix}.spike.sam || { echo -e "\n\033[7mSAMTOOLS SORT hap 1 failed! \033[0m";exit 1; } > ${hap1PreTumourBamPrefix}.spike.bam
${SAMTOOLS} index ${hap1PreTumourBamPrefix}.spike.bam || { echo -e "\n\033[7mSAMTOOLS INDEX hap 1 failed! \033[0m";exit 1; }

date | tr '\012' ':'
echo " Sort, index and convert to ${hap2PreTumourBamPrefix}.spike.bam."
${SAMTOOLS} sort ${hap2PreTumourBamPrefix}.spike.sam || { echo -e "\n\033[7mSAMTOOLS SORT hap 2 failed! \033[0m";exit 1; } > ${hap2PreTumourBamPrefix}.spike.bam
${SAMTOOLS} index ${hap2PreTumourBamPrefix}.spike.bam || { echo -e "\n\033[7mSAMTOOLS INDEX hap 2 failed! \033[0m";exit 1; }

# Remove redundant SAMs.
rm ${hap1PreTumourBamPrefix}.spike.sam ${hap2PreTumourBamPrefix}.spike.sam
