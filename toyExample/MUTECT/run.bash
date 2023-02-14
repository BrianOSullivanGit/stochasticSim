#!/bin/bash

# Use standard tool versions linked in with this toolkit.
# User may use other versions in their path on their local machine if they wish.
# They will probably work, but as I havent tested it, I can not guarantee it.
# The path to the set of tools that have been tested is specified in a file
# tool.path that is located in ..../stochasticSim-x.y/bin.

DT_STRING=`date -Iseconds`

TUMOUR_WITH_RG=`echo ${1} | sed -e "s/\.bam$/.rg.bam/1" -e "s/.*\///1"`
NORMAL_WITH_RG=`echo ${2} | sed -e "s/\.bam$/.rg.bam/1" -e "s/.*\///1"`

TUMOUR_BAM=`echo ${1} | sed "s/.*\///1"`
NORMAL_BAM=`echo ${2} | sed "s/.*\///1"`

VCF_NAME=`echo ${1} | sed -e "s/.*\///1" -e 's/^[TN]_//1' -e 's/\.bam$//1'`


date | tr '\012' ' '
echo "Adding read groups.."
samtools addreplacerg -r ID:HG00110 -r DT:${DT_STRING} -r SM:T_HG00110 ${1} -o ${TUMOUR_WITH_RG}
DT_STRING=`date -Iseconds`
samtools addreplacerg -r ID:HG00110 -r DT:${DT_STRING} -r SM:N_HG00110 ${2} -o ${NORMAL_WITH_RG}


LATEST_MUTECT_JAR="../gatk-4.2.2.0/gatk-package-4.2.2.0-local.jar"

samtools index ${TUMOUR_WITH_RG}
samtools index ${NORMAL_WITH_RG}

TUMOUR_SAMPLE=`${SAMTOOLS} view -H ${TUMOUR_WITH_RG} | grep -m1 '@RG.*SM:' | sed -e 's/@RG.*SM://1' -e 's/[[:space:]]*$//1'`
NORMAL_SAMPLE=`${SAMTOOLS} view -H ${NORMAL_WITH_RG} | grep -m1 '@RG.*SM:' | sed -e 's/@RG.*SM://1' -e 's/[[:space:]]*$//1'`

# Run mutect.
java -Xmx16g -jar ${LATEST_MUTECT_JAR} \
             Mutect2 \
            -R ../GRCh38.d1.vd1.chr19.fa \
            -I ${TUMOUR_WITH_RG} \
            -I ${NORMAL_WITH_RG} \
            -tumor ${TUMOUR_SAMPLE} \
            -normal ${NORMAL_SAMPLE} \
            -panel-of-normals ../1000g_pon.hg38.chr19.vcf.gz \
            -tumor-lod-to-emit 0 \
            -O ${VCF_NAME}.vcf.gz 

# Now filter calls (to annotate PASS variants etc).
java -Xmx16g -jar ${LATEST_MUTECT_JAR} \
      FilterMutectCalls \
      -R ../GRCh38.d1.vd1.chr19.fa \
      -V ${VCF_NAME}.vcf.gz \
      -O ${VCF_NAME}.filtered.vcf.gz


# Now create the map of called somatic mutations with the ground truth.
# This maps loci in the mutect VCF with their corresponding location in the donor personalised, phased genome
# and returns details about whether it is a real somatic mutation, its ground truth AF in the tumour etc.
# Map what gets called and what gets filtered to the ground truth.
../../bin/prepGtMap.bash \
                              ${VCF_NAME}.filtered.vcf.gz \
                              ../liftover_X1_HG00110.condense.txt \
                              ../liftover_X2_HG00110.condense.txt \
                              ${3} \
                              ${4} 2>&1 | tee summary.txt
echo "done."

