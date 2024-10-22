#!/bin/bash

# prepGtMap.bash /common/user1/SIM_PIPELINE/100x/MUTECT/HG00110.coding_exons_hg38_100x_76bp.realn.phased.filtered.vcf.gz \
#                /common/user1/1000GENOMES/TEST3/liftover_X1_HG00110.condense.txt \
#                /common/user1/1000GENOMES/TEST3/liftover_X2_HG00110.condense.txt \
#                /common/user1/SIM_PIPELINE/100x/T_X1_HG00110.coding_exons_hg38_100x_76bp.truth.vcf \
#                /common/user1/SIM_PIPELINE/100x/T_X2_HG00110.coding_exons_hg38_100x_76bp.truth.vcf
#                /common/user1/Reference/HG00110.vcf


#
# Use standard tool versions linked in with this toolkit.
# The user should have sourced the tool.path file to pick up the versions of tools
# on which this script depends. Otherwise the script will run with whatever versions it finds
# in the path. The user may already have other verions of these tools on their local machine
# and wish to run this script without sourcing the tool.path file.
# That will probably work, but as I havent tested it, I can not guarantee it.
# To use the set of tools against which this script has been tested please
# source the 'tool.path' file in the simulation framework bin directory before running it.
if [ -z ${TWOWAYLIFTOVER+xyz} ]; then TWOWAYLIFTOVER=`which 2wayLiftover`; fi
if [ -z ${GTMAPPER+xyz} ]; then GTMAPPER=`which gtMapper`; fi
if [ -z ${DATAMASH+xyz} ]; then DATAMASH=`which datamash`; fi
if [ -z ${BEDTOOLS+xyz} ]; then BEDTOOLS=`which bedtools`; fi

# Note:
# Do we need to reorder by haplotype, then reference coordinate order?
# Will these virtually always be the same?
# (decision, yes, so no action..) 

# Find target loci (loci of interest in the VCF for which we want fo pull out the ground truth) in haplotype 1.
gzip -cd ${1} | egrep -v '^#' | egrep '^chr[0-9]*[[:space:]]|^chr[XY][[:space:]]' | awk -v liftoverFile="${2}" '{if($4=="G" || $4=="C" || $4=="A" || $4=="T") if($5=="G" || $5=="C" ||$5=="A" || $5=="T") print "${TWOWAYLIFTOVER} r "$1":"$2" "liftoverFile}' > tmp.cmds.$$


date | tr '\012' ':'
echo "Lifting over target loci for haplotype 1, using ${2}"
# Remember at this point these loci are in reference coordinate order.
source tmp.cmds.$$ > targetsHap1.loci

# Cross check target loci against the GT VCF.
(${GTMAPPER} ${4} targetsHap1.loci || { echo -e "\n\033[7mGTMAPPER failed with haplotype 1 ! \033[0m" 1>&2; kill -HUP $$; }) | awk '{ gsub(".*AF=","",$8);if($4!=".") print $1":"$2"\t"$4">"$5"\t"$7"\t"$8; else print $1":"$2"\t.\t.\t." }' > gtMapper.hap1

# Repeat for Haplotype 2.

# Find target loci (loci of interest in the VCF for which we want fo pull out the ground truth) in haplotype 2.
gzip -cd ${1} | egrep -v '^#' | egrep '^chr[0-9]*[[:space:]]|^chr[XY][[:space:]]' | awk -v liftoverFile="${3}" '{if($4=="G" || $4=="C" || $4=="A" || $4=="T") if($5=="G" || $5=="C" ||$5=="A" || $5=="T") print "${TWOWAYLIFTOVER} r "$1":"$2" "liftoverFile}' > tmp.cmds.$$


date | tr '\012' ':'
echo "Lifting over target loci for haplotype 2, using ${3}"
# Remember at this point these loci are in reference coordinate order.
source tmp.cmds.$$ > targetsHap2.loci

# Cross check target loci against the GT VCF.
(${GTMAPPER} ${5} targetsHap2.loci || { echo -e "\n\033[7mGTMAPPER failed with haplotype 2 ! \033[0m" 1>&2; kill -HUP $$; }) | awk '{ gsub(".*AF=","",$8);if($4!=".") print $1":"$2"\t"$4">"$5"\t"$7"\t"$8; else print $1":"$2"\t.\t.\t." }' > gtMapper.hap2

# Get the filter fields and caller inferred allele frequencies from the VCF for all SBSs (again in ref. coordinate order).
# NB: This assumes the tumour format attributes are located in field 11 of the VCF record!!!
# (ie., this is a tumour normal pair from a GATK generated VCF with > v4.0.
# If it is not then your AF values man not be correct!!)
#
# 7/4/23, modified to allow this to work with Mutect2 tumout-only VCF files
gzip -cd ${1} | egrep -v '^#' | awk '{if($4=="G" || $4=="C" || $4=="A" || $4=="T") if($5=="G" || $5=="C" ||$5=="A" || $5=="T") {split($9,format,":");if(NF == 10) split($10,formatContentsTumour,":"); else split($11,formatContentsTumour,":"); for(i in format) {formatAttributesTumour[format[i]]=formatContentsTumour[i]}; print $1":"$2"\t"formatAttributesTumour["AF"]"\t"$7} }' > targets.filters.ref

# Finally combine the two into one file, removing redundant columns.
paste targets.filters.ref gtMapper.hap1 gtMapper.hap2 > gtMapper.hap.ref

# Now stats report.
hap1Burden=`grep -v '^#' ${4} | egrep -c '(PASS|MASK|UNDETECT|NO_COVER)'`
hap2Burden=`grep -v '^#' ${5} | egrep -c '(PASS|MASK|UNDETECT|NO_COVER)'`
burden=$((hap1Burden+hap2Burden))

hap1TotalUndetected=`grep -v '^#' ${4} | egrep -c '(UNDETECT|NO_COVER)'`
hap2TotalUndetected=`grep -v '^#' ${5} | egrep -c '(UNDETECT|NO_COVER)'`
totalUndetected=$((hap1TotalUndetected+hap2TotalUndetected))

hap1NoCoverage=`grep -v '^#' ${4} | egrep -c '(NO_COVER)'`
hap2NoCoverage=`grep -v '^#' ${5} | egrep -c '(NO_COVER)'`
noCoverage=$((hap1NoCoverage+hap2NoCoverage))


cat gtMapper.hap.ref |awk '{if($3=="PASS") print $0}' | cut -f 1,2,4-|egrep -v '(PASS|MASKED)' > fp.$$.list

egrep -v '^#' ${4} ${5}| egrep 'PASS|MASK' > delme.$$

totalNumSomaticsInSample=`wc -l delme.$$|cut -d' ' -f1`
totalNumMaskedInSample=`grep -c MASK delme.$$`

callerTotalPassed=`awk 'BEGIN{total=0}{if($3=="PASS") total++}END{print total}' gtMapper.hap.ref`

printf "
NGS:

The target contained %s true somatic variants.
%s of this number were not sequenced by the NGS process (ie., alt. allele depth = 0)
\n%s of those had zero coverage at variant locus (alt. allele depth = ref. allele depth = 0).

DNA containing alternative alleles from  %s somatic, single base substitutions (SBS) loci were picked up by NGS.

Of this %s sites were at least patrially masked by sequencing error.

Somatic Variant Calling:

The caller detected (annotated as \"PASS\") %s SBS in its output.

" "$burden" "$totalUndetected" "$noCoverage" "$totalNumSomaticsInSample" "$totalNumMaskedInSample" "$callerTotalPassed"



# If there are any false positives, produce a report.
if [ -s fp.$$.list ]
 then


   # Are we dealing with a tumour only VCF?
   numVcfFields=`zcat  ${1}| egrep -m1 '^#CHROM' | wc -w`

   # Only 10 fields in a VCF record means (in the context of somatic var. calling)  it's tumour only..
   if [ $numVcfFields -eq 10 ]; then
     # Tumour only.

     # There's no matched normal here so there will be no normal_artifact filter.
     # As a result there will be a tonne of germline false positives in the Mutect VCF output.
     # It is not informative to include them all in the summary output.
     # Just record the total and remove them from the list of false positives.
     gzip -cd ${1} | awk 'BEGIN{print "##fileformat=VCFv4.2"}{if($7=="PASS" && $4 ~ /^[GCAT]$/ && $5 ~ /^[GCAT]$/) print $0 }' > $$.tmp.vcf
     ${BEDTOOLS} intersect -a $$.tmp.vcf -b ${6} | awk '{print $1":"$2"\t"}' > n.a.tns.$$.lst
     germlineFalsePositives=`cat n.a.tns.$$.lst|wc -l`
     
     # Print out the totals..
     fpTotal=`wc -l fp.$$.list|cut -d' ' -f1`
     seqErrTotal=`grep SEQ_ERROR fp.$$.list | egrep -v -f n.a.tns.$$.lst| wc -l`
     alignmentIssue=`cat fp.$$.list|egrep -v -f n.a.tns.$$.lst|awk '{ if($4==".") if($5==".") if($6==".")if($8==".") if($9==".") if($10==".") print $0}' | wc -l`
     printf "\n\nFalse Positives.\n%s sites were incorrectly labeled as true somatic SBSs by the caller.\n%s of these were probably caused by incorrect base calls during the NGS process.\n%s of the remaining total were most likely caused by alignment issues post NGS.\n%s were germline artefacts mistakenly classified as true somatic variants.\n" "$fpTotal" "$seqErrTotal" "$alignmentIssue" "$germlineFalsePositives"
     printf "A list of false positives (artefacts incorrectly annotated as \"PASS\" by Mutect2) without the germline false positives, is shown below.\n"
     printf "NB: In the table, the true allele frequency in the haplotype (ie.,HAP[12]_true_AF) is twice
what is observed in the final (diploid) sequencing data.\nThis is because the final sequencing data is made up of sequencing data from both haplotypes combined together.\n\n"
     printf "VCF_locus\tVCF_AF\tHAP1_locus\tHAP1_SBS\tHAP1_ground_truth\tHAP1_true_AF\tHAP2_locus\tHAP2_SBS\tHAP2_ground_truth\tHAP2_true_AF\n"
     egrep -v -f n.a.tns.$$.lst  fp.$$.list

     rm $$.tmp.vcf
   else
     # Tumour and matched normal.
     # Remember, information about true negatives that are germline in origin and are filtered as
     # normal_artifact will not be found in the ground truth map file.
     # All germline information is recorded in the germline VCF file that was originally used to
     # create the personalised reference etc. (for example HG00110.vcf).
     # Most variants are obviously germline in origin are removed by Mutect2 without continuing through
     # the filtering process or recording them in the VCF to save processing time.
     # This significantly reduces the number of true negatives among variants listed as
     # normal_artifact. However we still need to check here just in case if there are any and remove them
     # from our normal_artifact false negative calculation.
     #
     # Get a list of germline origin normal_artifact true negatives and remove them from our false negatives list.
     # We will use this later on below when printing out the false negative list.
     gzip -cd ${1} | awk 'BEGIN{print "##fileformat=VCFv4.2"}{if($7=="normal_artifact" && $4 ~ /^[GCAT]$/ && $5 ~ /^[GCAT]$/) print $0 }' > $$.tmp.vcf
     ${BEDTOOLS} intersect -a $$.tmp.vcf -b ${6} | awk '{print $1":"$2"\t"}' > n.a.tns.$$.lst
     rm $$.tmp.vcf
     
     # Print out the totals..
     fpTotal=`wc -l fp.$$.list|cut -d' ' -f1`
     seqErrTotal=`grep SEQ_ERROR fp.$$.list | wc -l`
     alignmentIssue=`cat fp.$$.list|awk '{ if($4==".") if($5==".") if($6==".")if($8==".") if($9==".") if($10==".") print $0}' | wc -l`
     printf "\n\nFalse Positives.\n%s sites were incorrectly labeled as true somatic SBSs by the caller.\n%s of these were probably caused by incorrect base calls during the NGS process.\n%s of the remaining total were most likely caused by alignment issues post NGS.\n" "$fpTotal" "$seqErrTotal" "$alignmentIssue"
     printf "A list of false positives (artefacts incorrectly annotated as \"PASS\" by Mutect2) is shown below.\n"
     printf "NB: In the table, the true allele frequency in the haplotype (ie.,HAP[12]_true_AF) is twice
what is observed in the final (diploid) sequencing data.\n
This is because the final sequencing data is made up of sequencing data from both haplotypes combined together.\n\n"
     printf "VCF_locus\tVCF_AF\tHAP1_locus\tHAP1_SBS\tHAP1_ground_truth\tHAP1_true_AF\tHAP2_locus\tHAP2_SBS\tHAP2_ground_truth\tHAP2_true_AF\n"

     # Finally print out the list of false positives in the matched tumour normal VCF.
     cat fp.$$.list
   fi
else
   printf "\nThere were no false positives in the caller VCF output.\n"
fi

# Compile a report of filtered false negatives.
cat gtMapper.hap.ref|awk '{if($3!="PASS") print $0}' | egrep '(PASS|MASK)'|  sed 's/\t[^\t]*;[^\t]*\t/\tmultiple_filter_failures\t/1' | ${DATAMASH} --sort  --group 3 count 3 >  delme.$$

if [ -s delme.$$ ]; then
  totalFiltered=`cat delme.$$| ${DATAMASH} sum 2`
else
  totalFiltered="0"
fi
rm n.a.tns.$$.lst


printf "
False Negatives.

%s true somatic SBSs were not detected by the caller.
%s of these were incorrectly filetered as artifacts.
A breakdown of number of filtered false negative SBSs per filter type is shown in the associated plot.
The remaining true somatic SBSs were not listed in the caller VCF output.

" "$((burden + fpTotal - callerTotalPassed))" "$totalFiltered"

filtr=`cut -f1 delme.$$ | sed -e 's/^/"/1' -e 's/$/",/1' | tr -d '\012' | sed -e 's/^/c(/1' -e 's/,$/)/1'`
value=`cut -f2 delme.$$ | tr '\012' ',' | sed -e 's/^/c(/1' -e 's/,$/)/1'`

# Get output pdf name..
outputpdf=`echo ${1} | sed -e 's/.*\///1' -e 's/\.vcf.*/.fn.pdf/1'`

rscriptPath=`command -v Rscript` || { rscriptPath="/usr/bin/Rscript"; echo -e "\n\033[7mWARNING: RSCRIPT not found. False Neg. plot will not run on this machine. Consider installing R ! \033[0m" 1>&2; }


printf '#!%s
library(ggplot2)
library(ggrepel)
library(vctrs)
library(forcats)


# Define utility functions..
# This is a bit of a hack to get around install problems with tidyverse install on cluster.
# It provides a (considerably) more lightweight alternative.
tmpLead <- function (x, n = 1L, default = NA, order_by = NULL, ...) 
{
    if (!is.null(order_by)) {
        return(with_order(order_by, lead, x, n = n, default = default))
    }
    if (length(n) != 1 || !is.numeric(n) || n < 0) {
        bad_args("n", "must be a nonnegative integer scalar, ","not friendly_type_of(n) of length length(n).")
    }
    if (n == 0) 
        return(x)
    if (vec_size(default) != 1L) {
        abort(glue("default must be size 1, not size {vec_size(default)}"))
    }
    xlen <- vec_size(x)
    n <- pmin(n, xlen)
    inputs <- vec_cast_common(default = default, x = x)
    vec_c(vec_slice(inputs$x, -seq_len(n)), vec_rep(inputs$default, 
        n))
}

tmpCheck_factor <- function(f) { 
   if (is.character(f)) { 
     factor(f) 
   } else if (is.factor(f)) { 
     f 
   } else { 
     stop("`f` must be a factor (or character vector).", call. = FALSE) 
   } 
 } 


tmpFct_inorder <- function (f, ordered = NA) 
{
    f <- tmpCheck_factor(f)
    idx <- as.integer(f)[!duplicated(f)]
    idx <- idx[!is.na(idx)]
    lvls_reorder(f, idx, ordered = ordered)
}
value = %s           
csum = rev(cumsum(rev(value)))
pos = value/2 + tmpLead(csum, 1)
pos[which(is.na(pos))] = value[which(is.na(pos))]/2

df2 <- data.frame(
  group = %s,
  value = value,
  csum = csum,
  pos = pos
  ) 
           
pdf(file = "%s",   # where you want to save the file
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
ggplot(df2, aes(x = "" , y = value, fill = tmpFct_inorder(group))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_label_repel(data = df2,
                     aes(y = pos, label = value),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Filtered false negatives")) +
    theme_void()
dev.off()\n' "$rscriptPath" "$value" "$filtr" "${outputpdf}" > plotGtPieChart.R

# VAF plots
printf '#!%s

pdf(file = "%s",   # where you want to save the file
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches
    
# Shell commands to parse out the ground truth allele frequencies from spike-in files.
# If they are not found in the parent directory it ignores it and moves on.
groundTruthAlleleFrequenciesParseCommand="cat ../*.spike| awk %c{ print $4/2 }%c"
try({x = unlist(read.table(pipe(groundTruthAlleleFrequenciesParseCommand)));hist(x,xlim=c(0,1),breaks=100, main="Ground truth mutant allele frequencies", ylab="Number of mutations", xlab="Allele frequency")}, silent=TRUE)

# Shell commands to parse out allele frequencies from VCF output.
vcfAlleleFrequenciesParseCommand="gzip -cd *.filtered.vcf.gz | egrep -v %c^#%c | awk %c{if($7==\\"PASS\\" && $4 ~ /^[GCAT]$/ && $5 ~ /^[GCAT]$/) {split($9,format,\\":\\"); if(NF == 10) split($10,formatContentsTumour,\\":\\"); else split($11,formatContentsTumour,\\":\\"); for(i in format){formatAttributesTumour[format[i]]=formatContentsTumour[i]}; print formatAttributesTumour[\\"AF\\"]} }%c"

x = unlist(read.table(pipe(vcfAlleleFrequenciesParseCommand)))
hist(x,xlim=c(0,1),breaks=100, main="Mutant allele frequencies as estimted by Mutect2", ylab="Number of mutations", xlab="Allele frequency")

dev.off()

\n' "${rscriptPath}" "plotsVAF.pdf" "'" "'" "'" "'" "'" "'" > plotsVAF.R

# Remove temp files
rm targetsHap1.loci gtMapper.hap1 tmp.cmds.$$ targetsHap2.loci gtMapper.hap2 targets.filters.ref fp.$$.list delme.$$

# Create pie chart output.
chmod u+x plotGtPieChart.R plotsVAF.R
command -v Rscript > /dev/null && { ./plotGtPieChart.R; ./plotsVAF.R; }
exit 0
