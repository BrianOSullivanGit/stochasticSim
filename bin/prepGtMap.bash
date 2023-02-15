
# prepGtMap.bash /data/bosullivan/SIM_PIPELINE/100x/MUTECT/HG00110.coding_exons_hg38_100x_76bp.realn.phased.filtered.vcf.gz \
#              /data/bosullivan/1000GENOMES/TEST3/liftover_X1_HG00110.condense.txt \
#              /data/bosullivan/1000GENOMES/TEST3/liftover_X2_HG00110.condense.txt \
#              /data/bosullivan/SIM_PIPELINE/100x/T_X1_HG00110.coding_exons_hg38_100x_76bp.truth.vcf \
#              /data/bosullivan/SIM_PIPELINE/100x/T_X2_HG00110.coding_exons_hg38_100x_76bp.truth.vcf
#


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
if [ -z ${GTMAPPER+xyz} ]; then GTMAPPER=`which gtMapper"`; fi
if [ -z ${DATAMASH+xyz} ]; then DATAMASH=`which datamash"`; fi

# Note:
# Do we need to reorder by haplotype, then reference coordinate order?
# Will these virtually always be the same?
# (decision, yes, so no action..) 

# Find target loci (loci of interest in the VCF for which we want fo pull out the ground truth) in haplotype 1.
zcat ${1} | egrep -v '^#' | awk -v liftoverFile="${2}" '{if($4=="G" || $4=="C" || $4=="A" || $4=="T") if($5=="G" || $5=="C" ||$5=="A" || $5=="T") print "${TWOWAYLIFTOVER} r "$1":"$2" "liftoverFile}' > tmp.cmds.$$


date | tr '\012' ':'
echo "Lifting over target loci for haplotype 1, using ${2}"
# Remember at this point these loci are in reference coordinate order.
source tmp.cmds.$$ > targetsHap1.loci

# Cross check target loci against the GT VCF.
${GTMAPPER} ${4} targetsHap1.loci | awk '{ gsub(".*AF=","",$8);if($4!=".") print $1":"$2"\t"$4">"$5"\t"$7"\t"$8; else print $1":"$2"\t.\t.\t." }' > gtMapper.hap1


# Repeat for Haplotype 2.

# Find target loci (loci of interest in the VCF for which we want fo pull out the ground truth) in haplotype 2.
zcat ${1} | egrep -v '^#' | awk -v liftoverFile="${3}" '{if($4=="G" || $4=="C" || $4=="A" || $4=="T") if($5=="G" || $5=="C" ||$5=="A" || $5=="T") print "${TWOWAYLIFTOVER} r "$1":"$2" "liftoverFile}' > tmp.cmds.$$


date | tr '\012' ':'
echo "Lifting over target loci for haplotype 2, using ${3}"
# Remember at this point these loci are in reference coordinate order.
source tmp.cmds.$$ > targetsHap2.loci

# Cross check target loci against the GT VCF.
${GTMAPPER} ${5} targetsHap2.loci | awk '{ gsub(".*AF=","",$8);if($4!=".") print $1":"$2"\t"$4">"$5"\t"$7"\t"$8; else print $1":"$2"\t.\t.\t." }' > gtMapper.hap2

# Get the filter fields and caller inferred allele frequencies from the VCF for all SBSs (again in ref. coordinate order).
# NB: This assumes the tumour format attributes are located in field 11 of the VCF record!!!
# (ie., this is a tumour normal pair from a GATK generated VCF with > v4.0.
# If it is not then your AF values man not be correct!!)
zcat ${1} | egrep -v '^#' | awk '{if($4=="G" || $4=="C" || $4=="A" || $4=="T") if($5=="G" || $5=="C" ||$5=="A" || $5=="T") {split($9,format,":");split($11,formatContentsTumour,":"); for(i in format) {formatAttributesTumour[format[i]]=formatContentsTumour[i]}; print $1":"$2"\t"formatAttributesTumour["AF"]"\t"$7} }' > targets.filters.ref

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

fpTotal=`wc -l fp.$$.list|cut -d' ' -f1`

seqErrTotal=`grep SEQ_ERROR fp.$$.list | wc -l`

alignmentIssue=`cat fp.$$.list|awk '{ if($4==".") if($5==".") if($6==".")if($8==".") if($9==".") if($10==".") print $0}' | wc -l`

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

False Positives.
%s sites were incorrectly labeled as true somatic SBSs by the caller.
%s of these originated from incorrect base calls during the NGS process.
%s originated from likely alignment issues post NGS.

" "$burden" "$totalUndetected" "$noCoverage" "$totalNumSomaticsInSample" "$totalNumMaskedInSample" "$callerTotalPassed" "$fpTotal" "$seqErrTotal" "$alignmentIssue"

# If there are any false positives, produce a report.
if [ -s fp.$$.list ]
 then
   printf "A list of false positives (artefacts incorrectly annotated as \"PASS\" by Mutect2) is shown below.\n"
   printf "NB: In the table, the true allele frequency in the haplotype (ie.,HAP[12]_true_AF) is twice
what is observed in the final (diploid) sequencing data.
This is because the final sequencing data is made up of sequencing data from both haplotypes combined together.\n\n"
   printf "VCF_locus\tVCF_AF\tHAP1_locus\tHAP1_SBS\tHAP1_ground_truth\tHAP1_true_AF\tHAP2_locus\tHAP2_SBS\tHAP2_ground_truth\tHAP2_true_AF\n"

   cat fp.$$.list
fi


cat gtMapper.hap.ref|awk '{if($3!="PASS") print $0}' | egrep '(PASS|MASK)'| sed 's/\t[^\t]*;[^\t]*\t/\tmultiple_filter_failures\t/1' | ${DATAMASH} --sort  --group 3 count 3 >  delme.$$

totalFiltered=`cat delme.$$| ${DATAMASH} sum 2`

printf "
False Negatives.

%s true somatic SBSs picked up by NGS were not detected by the caller.
%s of these were incorrectly filetered as artifacts.
A breakdown of number of filtered false negative SBSs per filter type is shown in the associated plot.
The remainding true somatic SBSs were not listed in the caller VCF output.

" "$((burden - callerTotalPassed))" "$totalFiltered"

filtr=`cut -f1 delme.$$ | sed -e 's/^/"/1' -e 's/$/",/1' | tr -d '\012' | sed -e 's/^/c(/1' -e 's/,$/)/1'`
value=`cut -f2 delme.$$ | tr '\012' ',' | sed -e 's/^/c(/1' -e 's/,$/)/1'`

# Get output pdf name..
outputpdf=`echo ${1} | sed -e 's/.*\///1' -e 's/\.vcf.*/.fn.pdf/1'`


printf '#!/usr/bin/Rscript
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
    width = 15, # The width of the plot in inches
    height = 5) # The height of the plot in inches
ggplot(df2, aes(x = "" , y = value, fill = tmpFct_inorder(group))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_label_repel(data = df2,
                     aes(y = pos, label = value),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Filtered false negatives")) +
    theme_void()
dev.off()\n' "$value" "$filtr" "${outputpdf}"> plotGtPieChart.R

# Remove temp files
rm targetsHap1.loci gtMapper.hap1 tmp.cmds.$$ targetsHap2.loci gtMapper.hap2 targets.filters.ref fp.$$.list delme.$$

# Create pie chart output.
chmod u+x plotGtPieChart.R
./plotGtPieChart.R
