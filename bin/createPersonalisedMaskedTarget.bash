#!/bin/bash
# STEP 1 Create personalised, masked target.
# Outputs a pair of phased personalised reference fastas plus a pair of liftover files
# to map between them and the standard reference supplied that was supplied as input.
#
# createPersonalisedMaskedTarget.bash <ID prefix> <sex, m|f> <standard reference genome> <phased germline VCF> <target BED> <flanking region, if required>

# ../createPersonalisedMaskedTarget.bash HG00110 F /data/bosullivan/TMB/REF/GRCh38.d1.vd1.fa /data/bosullivan/1000GENOMES/TEST7/HG00110.vcf coding_exons_hg38.bed 100

# Use standard tool versions linked in with this toolkit.
# User may use other versions in their path on their local machine if they wish.
# They will probably work, but as I havent tested it, I can not guarantee it.
# The path to the set of tools that have been tested is specified in a file
# tool.path that is located in ..../stochasticSim-x.y/bin.

if [ -z ${BEDTOOLS+xyz} ]; then BEDTOOLS=`which bedtools`; fi
if [ -z ${CREATEDONORGENOME+xyz} ]; then CREATEDONORGENOME=`which createDonorGenome`; fi
if [ -z ${LIFTOVER+xyz} ]; then LIFTOVER=`which liftover`; fi
if [ -z ${TARGETREF+xyz} ]; then TARGETREF=`which targetRef`; fi

# If no flanking supplies, set it to zero so it will be ignored.
if [ -z "$6" ]
  then
    flankingLength=0
else
    flankingLength=${6}
fi


refPrefix=`echo ${3} | sed -e 's/.*\///1' -e 's/\.fa$//1'`
echo $refPrefix


bedPrefix=`echo ${5} | sed -e 's/.*\///1' -e 's/\.bed.*$//1'`
echo $bedPrefix


# Subset standard reference fasta passed in here, keeping only chr1-22XY (ie., everything before chrM).
date | tr '\012' ':'
cat /dev/null > ${refPrefix}.chr1-22XY.fa
echo ' Subset reference to contain only chr1-22XY.'
(seq 1 1 22|sed 's/^/chr/1';echo chrX;echo chrY) | while read line; do samtools faidx ${3} ${line} 2>/dev/null > $$.fa && cat $$.fa >> ${refPrefix}.chr1-22XY.fa; done
rm $$.fa


# Sort the bed (in same coordinate order as the corresponding reference file), keeping only chr1-22XY
cat ${5} |egrep '^chr[1-9][[:space:]]|^chr[1-2][0-9][[:space:]]' | sort -k1.4,1.5n -k2,2n -k3,3n|cut -f1-3 > ${bedPrefix}.sorted.bed
cat ${5} |egrep '^chrX[[:space:]]' | sort -k2,2n -k3,3n|cut -f1-3 >> ${bedPrefix}.sorted.bed
cat ${5} |egrep '^chrY[[:space:]]' | sort -k2,2n -k3,3n|cut -f1-3 >> ${bedPrefix}.sorted.bed

# Merge overlapping ranges.
${BEDTOOLS} merge -i ${bedPrefix}.sorted.bed  > tmp.$$.bed

# Allow for off-target capture in flanking regions.
# TODO!!! Make sure the flanking region arg was set before doing this...
# TODO!!!! improve estimate here (most capture panels are about 80% on-target, 20% off).
# For now just tag on ${flankingLength}bp either side of the genomic range in specified in the BED record.

if [ ${flankingLength} -ne 0 ]
then
  cat tmp.$$.bed | awk -v flank=${flankingLength} '{ print $1"\t"($2-flank)"\t"($3+flank) }'  > tmp.$$.${flankingLength}bpFlanking.bed
else
  mv tmp.$$.bed tmp.$$.${flankingLength}bpFlanking.bed
fi

# Re-sort.
cat ./tmp.$$.${flankingLength}bpFlanking.bed|egrep '^chr[1-9][[:space:]]|^chr[1-2][0-9][[:space:]]' | sort -k1.4,1.5n -k2,2n -k3,3n|cut -f1-3 > tmp.$$.bed
cat ./tmp.$$.${flankingLength}bpFlanking.bed|egrep '^chrX[[:space:]]' | sort -k2,2n -k3,3n|cut -f1-3 >> tmp.$$.bed
cat ./tmp.$$.${flankingLength}bpFlanking.bed|egrep '^chrY[[:space:]]' | sort -k2,2n -k3,3n|cut -f1-3 >> tmp.$$.bed


# Re-merge overlapping ranges.
${BEDTOOLS} merge -i tmp.$$.bed | gzip > ${bedPrefix}.sorted.filtered.${flankingLength}bpFlanking.merged.bed.gz

# Clean as you go..
rm tmp.$$.bed tmp.$$.${flankingLength}bpFlanking.bed ${bedPrefix}.sorted.bed

# Generate personalised genome.
# Female ${1} for example.
# Note the corresponding VCF must already be in this directory
# (if not extract it with extractDonorVcfOneKG.sh)
# TODO!!!! Make 1000G subject configurable..

# X1 set
date | tr '\012' ':'
echo ' Create personalised donor genome. Reference X1 (haplotype 1) X1_<ID prefix>.fa.'
if [ "${2}" = "F" -o "${2}" = "f" ]
then
    rm X1_${1}.fa.fai 2>/dev/null; ${CREATEDONORGENOME} ${refPrefix}.chr1-22XY.fa ${4} "0|1" chrY > X1_${1}.fa;mv liftover.txt liftover_X1_${1}.txt;date
else
    # TODO!!!! Male, don't skip Y!!!!!(or skip Y here but not below????)
    echo "Not fully implemented yet!!!!!!!"
    rm X1_${1}.fa.fai 2>/dev/null; ${CREATEDONORGENOME} ${refPrefix}.chr1-22XY.fa ${4} "0|1" > X1_${1}.fa;mv liftover.txt liftover_X1_${1}.txt;date
fi


# X2 set / Y2 set if male..
date | tr '\012' ':'
echo ' Create personalised donor genome. Reference X2 (haplotype 2) X2_<ID prefix>.fa.'


if [ "${2}" = "F" -o "${2}" = "f" ]
then
    rm X2_${1}.fa.fai 2>/dev/null; ${CREATEDONORGENOME} ${refPrefix}.chr1-22XY.fa ${4} "1|0" chrY > X2_${1}.fa;mv liftover.txt liftover_X2_${1}.txt;date
else
    echo "Not fully implemented yet!!!!!!!"
    rm X2_${1}.fa.fai 2>/dev/null; ${CREATEDONORGENOME} ${refPrefix}.chr1-22XY.fa ${4} "1|0" > X2_${1}.fa;mv liftover.txt liftover_X2_${1}.txt;date
fi


# Create new target BEDs (X1 set & X2 set) to match personalised diploid genome.
date | tr '\012' ':'
echo ' Liftover target BED to align with personalised donor genome.'
if [ "${2}" = "F" -o "${2}" = "f" ]
then
    ${LIFTOVER}  liftover_X1_${1}.txt <(zcat ${bedPrefix}.sorted.filtered.${flankingLength}bpFlanking.merged.bed.gz) chrY | gzip > X1_${1}_${bedPrefix}.merged.bed.gz
    ${LIFTOVER}  liftover_X2_${1}.txt <(zcat ${bedPrefix}.sorted.filtered.${flankingLength}bpFlanking.merged.bed.gz) chrY | gzip > X2_${1}_${bedPrefix}.merged.bed.gz
else
    echo "Not fully implemented yet!!!!!!!"
    ${LIFTOVER}  liftover_X1_${1}.txt <(zcat ${bedPrefix}.sorted.filtered.${flankingLength}bpFlanking.merged.bed.gz) | gzip > X1_${1}_${bedPrefix}.merged.bed.gz
    ${LIFTOVER}  liftover_X2_${1}.txt <(zcat ${bedPrefix}.sorted.filtered.${flankingLength}bpFlanking.merged.bed.gz) | gzip > X2_${1}_${bedPrefix}.merged.bed.gz
fi


# Mask personalised genome target BEDs
# (ie., overwrite non matching regions with "N"'s so no reads are generated from those regions).
date | tr '\012' ':'
echo " Subset personalised genome to only include required targets.."
${TARGETREF} X1_${1}.fa <(zcat X1_${1}_${bedPrefix}.merged.bed.gz) > X1_${1}.${bedPrefix}.fa
${TARGETREF} X2_${1}.fa <(zcat X2_${1}_${bedPrefix}.merged.bed.gz) > X2_${1}.${bedPrefix}.fa
