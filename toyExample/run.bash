#!/bin/bash

target=`echo ${2} | sed 's/\.bed$//1'`
hapDepth=`bc -l <<< "scale=0;${1}/2"`

echo -e "\033[7mStart Simulaion for ${target} at ${1X} for donor HG00110\033[0m\n"

printf "Create a phased, personalised reference genome that is modified to reflect target region of interest.
All locations outside the target area will have their contents replaced with 'N' bases to avoid reads
being simulated from these regions. This command creates the pair of phased fasta reference files for donor
HG00110 together with other files for use later in the simulation.\n\n"
date | tr '\012' ':'
echo -e ">>>>>> Start \033[1m../bin/createPersonalisedMaskedTarget.bash HG00110 \\ \n                                                                                     F \\ \n                                                                                     GRCh38.d1.vd1.chr19.fa \\ \n                                                                                     HG00110.chr19.vcf \\ \n                                                                                     ${2}\033[0m"

../bin/createPersonalisedMaskedTarget.bash HG00110 F GRCh38.d1.vd1.chr19.fa HG00110.chr19.vcf ${2}

date | tr '\012' ':'
echo -e ">>>>>> Done. \033[1m../bin/createPersonalisedMaskedTarget.bash\033[0m"
echo

printf "Cleanup liftover that are used to map coordinates from standard to personalised diploid reference.\n\n"
date | tr '\012' ':'
echo -e ">>>>>> Start \033[1mCompress liftover files\033[0m" 

../bin/condenseLift liftover_X1_HG00110.txt > liftover_X1_HG00110.condense.txt
../bin/condenseLift liftover_X2_HG00110.txt > liftover_X2_HG00110.condense.txt
rm liftover_X1_HG00110.txt liftover_X2_HG00110.txt

date | tr '\012' ':'
echo -e ">>>>>> Done. \033[1mCompress liftover files\033[0m"
echo

printf "Create the spike-in config file.
This matches the frequency distribution of interest with a set of target loci for the spike-in config.
Target loci are chosen according to your research requirements.
They may include particular drivers of interest or variants within a specific set of
trinucleotide contexts, or simply chosen at random (as in this example).\n\n"
date | tr '\012' ':'
echo -e ">>>>>> Start \033[1mspike-in config files\033[0m"

# Match your frequency distribution with a set of target loci for the spike-in config file.
# Target loci are chosen according to your research requirements.
# They may include particular drivers of interest or variants within a specific set of
# trinucleotide contexts, or simply chosen at random (as in this example).
cat ${target}.bed | awk '{ i=$2; while(i!=$3) {print $1"\t"i+1"\t.";i++} }' |shuf -n  `wc -l ${target}.h1.854mutMB.freqs|sed 's/ .*//1'` > ${target}.h1.854mutMB.loci
cat ${target}.bed | awk '{ i=$2; while(i!=$3) {print $1"\t"i+1"\t.";i++} }' |shuf -n  `wc -l ${target}.h2.854mutMB.freqs|sed 's/ .*//1'` > ${target}.h2.854mutMB.loci

# Paste the two together to create the final spike-in config.
# Remember it is essencial that all VCFs, bed files and spike-in config files are sorted
# in chromosome and coordinate order before using them with the simulation framework.
paste ${target}.h1.854mutMB.loci ${target}.h1.854mutMB.freqs | sort -k1,1V -k2,2n > ${target}.h1.854mutMB.spike
paste ${target}.h2.854mutMB.loci ${target}.h2.854mutMB.freqs | sort -k1,1V -k2,2n > ${target}.h2.854mutMB.spike

date | tr '\012' ':'
echo -e ">>>>>> Done. \033[1mspike-in config\033[0m"
echo

printf "Simulate a matched normal and pre-tumour pair of BAM files.
Use the personalised target reference files (for HG00110) to create a phased, matched normal and pre-tumour
set of BAMs at the required depth of coverage.\n\n"

date | tr '\012' ':'
echo -e ">>>>>> Start \033[1m../bin/generatePhasedBams.bash X1_HG00110.${target}.fa \\ \n                                                                         X2_HG00110.${target}.fa \\ \n                                                                         76 \\ \n                                                                         \"-f ${hapDepth}\" \\ \n                                                                         180 \\ \n                                                                         ${hapDepth}x_76bp\033[0m"

../bin/generatePhasedBams.bash X1_HG00110.${target}.fa X2_HG00110.${target}.fa 76 "-f ${hapDepth}" 180 "${hapDepth}x_76bp"

date | tr '\012' ':'
echo -e ">>>>>> Done. \033[1m../bin/generatePhasedBams.bash\033[0m"

printf "\nSpike-in the required burden to the phased pre-tumour BAM files.\n\n"
date | tr '\012' ':'
echo -e ">>>>>> Start \033[1m../bin/spikeIn.bash T_X1_HG00110.${target}_${hapDepth}x_76bp.bam \\ \n                                                              T_X2_HG00110.${target}_${hapDepth}x_76bp.bam \\ \n                                                              X1_HG00110.${target}.fa \\ \n                                                              X2_HG00110.${target}.fa \\ \n                                                              ${target}.h1.854mutMB.spike \\ \n                                                              ${target}.h2.854mutMB.spike\033[0m"
../bin/spikeIn.bash T_X1_HG00110.${target}_${hapDepth}x_76bp.bam T_X2_HG00110.${target}_${hapDepth}x_76bp.bam X1_HG00110.${target}.fa X2_HG00110.${target}.fa ${target}.h1.854mutMB.spike  ${target}.h2.854mutMB.spike
date | tr '\012' ':'
echo -e ">>>>>> Done, \033[1mspike-in.\033[0m"

printf "\nRe-align against the standard reference and merge tumour and normal phased BAM sets to
create a tumour and normal pair for somatic variant calling.\n\n"

date | tr '\012' ':'
echo -e ">>>>>> Start \033[1m../bin/realignAndMerge.bash T_X1_HG00110.${target}_${hapDepth}x_76bp.spike.bam  \\ \n                                                                      T_X2_HG00110.${target}_${hapDepth}x_76bp.spike.bam \\ \n                                                                      GRCh38.d1.vd1.chr19.fa \\ \n                                                                      T_HG00110.${target}_${1}x_76bp\033[0m"
../bin/realignAndMerge.bash T_X1_HG00110.${target}_${hapDepth}x_76bp.spike.bam  T_X2_HG00110.${target}_${hapDepth}x_76bp.spike.bam GRCh38.d1.vd1.chr19.fa T_HG00110.${target}_${1}x_76bp
date | tr '\012' ':'
echo -e ">>>>>> Done, \033[1mrealign & merge, tumour BAMs.\033[0m"

date | tr '\012' ':'
echo -e ">>>>>> Start \033[1m../bin/realignAndMerge.bash N_X1_HG00110.${target}_${hapDepth}x_76bp.bam  \\ \n                                                                      N_X2_HG00110.${target}_${hapDepth}x_76bp.bam \\ \n                                                                      GRCh38.d1.vd1.chr19.fa \\ \n                                                                      N_HG00110.${target}_${1}x_76bp\033[0m"
../bin/realignAndMerge.bash N_X1_HG00110.${target}_${hapDepth}x_76bp.bam  N_X2_HG00110.${target}_${hapDepth}x_76bp.bam GRCh38.d1.vd1.chr19.fa N_HG00110.${target}_${1}x_76bp
date | tr '\012' ':'
echo -e ">>>>>> Done, \033[1mrealign & merge, normal BAMs.\033[0m"

printf "\nPerform somatic variant calling (Mutect2). Map the output against the ground truth.\n\n"

date | tr '\012' ':'
echo -e ">>>>>> Start \033[1mSomatic variant caller and map its output to the ground truth.\033[0m"
cd MUTECT
./run.bash ../T_HG00110.${target}_${1}x_76bp.realn.phased.bam ../N_HG00110.${target}_${1}x_76bp.realn.phased.bam ../T_X1_HG00110.${target}_${hapDepth}x_76bp.truth.vcf ../T_X2_HG00110.${target}_${hapDepth}x_76bp.truth.vcf 
date | tr '\012' ':'
echo -e ">>>>>> Done, \033[1msomatic variant caller and ground truth map.\033[0m"
echo
echo -e "\033[7mSimulaion for ${target} at ${1X}\033, donor HG00110 complete\033[0m"
echo
echo "Look in "${PWD}" for the simulations results."
echo "Along with the variant caller output, this directory will contain the following files"
echo
echo -e "\033[1mgtMapper.hap.ref:\033[0m     A tab seperated table mapping all entries in the caller filtered VCF output to their\n                      ground truth values in each haplotype."
echo -e "\033[1mplotGtPieChart.R/pdf:\033[0m Source and associated pdf output to plot a pie chart of all caller filtered false negatives."
echo -e "\033[1msummary.txt:\033[0m          An overall breakdown of where caller false positives/negatives occurred in this simulation\n                      and why."
echo
cd ../
echo -e "For further details on how to run this simulation and others like it look at the contents of \n"${PWD}"/run.bash"



