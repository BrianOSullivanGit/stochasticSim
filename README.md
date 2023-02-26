# Stochastic simulation framework<!-- omit in toc -->

*A set of software tools providing comprehensive and realistic simulation of tumour genomic sequencing data.*

This framework accounts for the random nature of genomic sequencing to accurately reflect the frequency profile captured in real tumour genomic sequencing data. The ground truth associated with the simulation output is recorded for every non-reference locus in the data. Every record in the VCF output by the caller (either PASS or filtered) is mapped back to its true location in the phased personalised donor record to definitively identify false positives / negatives and their causes.

Contact BrianOSullivan@yahoo.com with questions.

## Table of contents<!-- omit in toc -->
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Getting started](#getting-started)
- [Understanding simulation output](#understanding-simulation-output)
- [List of stochasticSim commands and their format](#list-of-stochasticsim-commands-and-their-formats)

## System requirements

This simulation framework runs on a 64bit, Mac OS X or GNU/Linux. If this does not cover your system it will not work. If you need a version for a different platform let me know (BrianOSullivan@yahoo.com). If you want to run Mutect2 as your somatic caller (as in toyExample included) you will also need Java 8 / JRE or SDK 1.8 installed (check your version with `java -version`). The `prepGtMap.bash` script included which plots the contribution of each GATK filter to the number of false negatives in the VCF output also requires `R` (script, `/usr/bin/Rscript`) together with libraries `ggplot2`, `ggrepel`, `vctrs` and `forcats` (see https://cran.r-project.org if you do not already have R installed).

The total disc space taken up by the framework, once all associated tools (samtools, bwa etc.) are installed locally and the toy example simulation run is about 1.7G. Approximately 0.7G of this is taken up by BAM files created during simulation. These files may be deleted once the somatic variant caller has run (unless you need to look into alignment artefacts/issues in caller output, in that case you will need them).

This install also builds the following publically available tools for use with this framework.

    bwa >0.7.17.
    samtools 1.13
    bedtools 2.29.2

They are installed locally, in the same directory in which `stochasticSim` was cloned and neither the install or simulation process require root access to run. The installation process will not modify your path or cause any conflict with any previous installation of these tools you may have on your system. You may if you wish use other versions of these tools with `stochasticSim` if you have them installed. They will probably work, however as I have not tested it, I can not guarantee it. Links to all scripts and binaries used in this simulation framework are included under the top level bin directory (`stochasticSimFramework/bin`). You can copy them into another directory in your path and run them from there if you prefer. The install script uses curl to download these dependencies. You must have curl installed on your system or the install will not work (ie., for ubuntu linux, `apt install curl`). There are issues when building bwa aligner version 0.7.17 with [some gcc versions (> 10)](https://github.com/lh3/bwa/pull/385) so we use a clone of the bwa main branch in the install which contains a fix.

To build `stochasticSim` and associated tool dependencies you will need a basic set of developer resources set up such as GNU make, C / C++ compiler, autoconf e.t.c. On many installations these will already be in place (particularly if they have been used to build tools previously). If not however you will need to install them. On Ubuntu Linux this may require the [build-essential package](https://packages.ubuntu.com/focal/build-essential) `sudo apt-get install build-essential`, on MacOS select "command line developer tools" from `xcode-select -install` (see the [Ubuntu Package archive](https://packages.ubuntu.com) or [Apple's Xcode site](https://developer.apple.com/xcode) for relevant information and downloads). These installations will require root access. Samtools (& HTSlib) also depend on a number of libraries. Again, these may already be on you system. On Ubuntu linux, for example, you can check with,

```
$ apt -qq list autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses-dev

autoconf/jammy,jammy,now 2.71-2 all [installed]
automake/jammy,jammy,now 1:1.16.5-1.3 all [installed]
gcc/jammy,now 4:11.2.0-1ubuntu1 amd64 [installed]
gcc/jammy 4:11.2.0-1ubuntu1 i386
libbz2-dev/jammy,now 1.0.8-5build1 amd64 [installed]
                :
                :
                :
```
and if necessary, install what you need with,
```
$ sudo apt-get install libncurses-dev
                :
              e.t.c
                :
```
If you run into an issue with "python not found" when building bedtools on ubuntu, consider installing python-is-python3 (`sudo apt-get install python-is-python3`, again it is likely the fix will already be in place). Finally, the install script will also download ART read simulator (mountrainier2016) binaries and the `gatk-package-4.2.2.0` jar file (containing MUTECT2 somatic variant caller, requires Java 8 / JRE or SDK 1.8).



## Installation

Before running this script, run the following commands.
Make a directory into which the set of tools which make up the simulation framework will be installed. Enter that directory
```
cd <your download target directory>
mkdir stochasticSimFramework
cd stochasticSimFramework
```

Now download and unpack the stochasticSim repository, either latest release or main branch from gitHub. For example using main,
```
unzip <download path>/stochasticSim-main.zip
```

Enter the `stochasticSim` directory. Run the install script, checking for any errors in the output and if required rectify any issues that may have occurred during installation.
```
cd stochasticSim
more ./install.bash
./install.bash
```

## Getting started
Once the install has completed successfully you can test it out on the toy example included.
First source the `tool.path` file to ensure the simulation uses
the required versions of associated software tools.
```
source <install path>/stochasticSimFramework/stochasticSim-main/bin/tool.path
```
Now, enter the toy example directory and run the simulation,
```
cd <install path>/stochasticSimFramework/stochasticSim-main/toyExample
./run.bash 50 chr19_500KB.bed
```
Depending on your hardware this will take about 8 minutes to run.

The sequencing data used with this framework is simulated from phased, personalised reference files containing all germline SNV and indels variation from a real donor (HG00110, from the 1000 genomes project, in the toy example included). This is achieved by creating two personalised reference fasta files, each notionally corresponding to the maternal and paternal set of chromosomes, based on the donor's phased germline VCF. These fasta files are then used to simulate two BAM files using the ART read simulator (configured with default empirical error profile and corresponding to 50x depth of coverage in total, for the toy example). In the case of the simulated tumour sample, these BAMs are used as a base to spike in a set of somatic variants at the required loci and true allele frequencies. The BAMs are then realigned against the standard reference and merged before input to the somatic variant calling pipeline of interest. At each stage in the process, information is recorded about the source and location of every non reference site as it arises in the data. This information is critical in evaluating variant caller output at the end of this process. The files that store this information and the format used is outlined in the next section.


## Understanding simulation output
Look at `<install path>/stochasticSimFramework/stochasticSim-main/toyExample/MUTECT` for the simulations results.
Along with the variant caller output, this directory will contain the following files

| Filename | Description |
| --- | --- |
| `gtMapper.hap.ref` | A tab separated table mapping all entries in the caller filtered VCF output to their ground truth values in each haplotype. |
| `plotGtPieChart.R/pdf` | Source and associated pdf output to plot a pie chart of all caller filtered false negatives. |
| `summary.txt` | An overall breakdown of where caller false positives/negatives occurred in this simulation and why. |
| `T_X1_<id>.truth.vcf` | A VCF recording the true origin of every site in the first phased BAM that does not match the personalised reference. |
| `T_X\|Y2_<id>.truth.vcf` | A VCF recording the true origin of every site in the second phased BAM that does not match the personalised reference. |

### Ground truth VCFs

The two ground truth VCFs files (one for each haplotype, *.truth.vcf) are created during somatic variant spike-in on BAM data that is perfectly aligned aginst the appropriate personalised haploid reference of the phased pair. The purpose of these files is to record the true reason behind all non-reference loci in the data so we can explain why somatic variants are missed or incorrectly called later in the pipeline when the variant calling has been run. The filter field, indicating the true reason for the non-reference site is enumerated as follows,

| Filter field | Description |
| --- | --- |
| SEQ_ERROR | The locus does not contain an alternative allele. The non reference base(s) is a result of sequence error during simulation |
| PASS | The locus contains a somatic variant and all reads supporting that variant were successfully spiked into the data |
| MASKED | The locus contains a somatic variant however one or more reads supporting that somatic variant were masked by sequence error |
| MASKED_OVL | The locus contains a somatic variant however one or more reads supporting that somatic variant were masked by sequence error in a read pair that overlapped at the target locus. |
| NO_COVERAGE | The locus contains a somatic variant however no coverage was obtained at this locus during simulation (DP=0)|
| UNDETECTED | The locus contains a somatic variant however no reads containing the alternative allele were detected during the simulation (AD for the spiked-in allele=0) |

For PASS and MASKED entries, the AF attribute of the INFO field will contain the true allele frequency of the spiked-in somatic variant in that haplotype. Note that, as both tumour haplotype BAM files are combined together after realignment to create the final tumour BAM, the allele frequency listed ground truth VCF for the haplotype in question is approximately twice its actual value in the final tumour data. The user should be aware that stochasticSim accounts for the randomness in the number of reads that contain the non-reference allele at somatic mutation sites. This means that in our simulations (as in real sequencing data) the true allele frequency and the allelic ratio (AD/DP) will not necessarily match (though, in general, particularly at reasonably high depth, they will be close). PASS and MASKED entries will always list the allele depth (AD) of the spiked-in allele second (directly after the reference allele depth) in the SAMPLE field. For SEQ_ERROR entries the AF entry in the INFO field will contain the allelic ratio of the first error "allele" (the one with the highest alt. allele depth). Similarly (for SEQ_ERROR) the depths of any other error alleles will be listed in descending order after the reference allele depth in the SAMPLE field. A sequencer error that (by chance) creates the same allele as the somatic variant at a MASKED locus will be listed separately in the SAMPLE field. For example the MASKED variant record from the selection of ground truth VCF entries below indicates that 3 reads containing the somatic 'A' allele were picked up by the sequencer. It also indicates that a sequencer error occurred in one or more other reads which also contained the somatic 'A' allele, causing it to be incorrectly read as a 'T'. Finally, another sequencer error occurred at the same locus causing one of the reads containing the reference 'C' allele to be incorrectly read as an 'A'. The vast majority of entries in the ground truth VCF files however will record non-reference sites caused by single read sequence error (SEQ_ERROR, AD=19,1 for example).


```
#CHROM     POS             ID      REF     ALT     QUAL    FILTER      INFO               FORMAT  SAMPLE
chr1       69073           .       A       G       .       SEQ_ERROR   DP=35;AF=0.02857   AD      34,1
chr1       1180426         .       A       T       .       PASS        DP=49;AF=0.05906   AD      46,3
chr1       3499182         .       A       C       .       UNDETECTED  DP=16;AF=0.01580   AD      16,0
chr1       111311810       .       C       A,T,A   .       MASKED      DP=39;AF=0.05643   AD      32,3,3,1
chr3       193301639       .       .       .       .       NO_COVERAGE     .       .
```

The ground truth VCFs are created prior to somatic variant calling, when the simulated data is perfectly aligned against the phased personalised donor reference. Issues that arise after the data is realigned against the standard reference and used as input to a somatic variant caller are recorded in the ground truth map file (`gtMapper.hap.ref`). This file, which maps each entry in the somatic caller VCF to its corresponding entries in each of the ground truth VCF files, is described in the following section.

### Ground truth map file

This file (`gtMapper.hap.ref`) maps each record in the (filtered VCF) somatic variant caller output to its corresponding ground truth in each haplotype (recorded in the ground truth VCFs). Fields 1 to 3 of each record in the ground truth map show the genomic coordinates (chromosome:locus, 1 base), allele frequency (as estimated by the caller) and FILTER field of the variant as it appears in the caller output VCF. The coordinates indicate its aligned location in the standard reference genome. Fields 4 to 7 show the corresponding genomic coordinates in the first haplotype, the single base substitution in the first haplotype (if any), FILTER field contents and allele frequency in the haplotype taken from the ground truth VCF for haplotype 1 (`T_X1_<id>.truth.vcf`). Recall that, as both haplotype BAM files are combined together after realignment to create the final tumour BAM, the allele frequency listed in the haplotype is approximately twice its actual value in the final tumour data. Finally, fields 8 to 11 show the corresponding genomic coordinates in the second haplotype, the single base substitution (if any), FILTER field contents and allele frequency in the haplotype taken from the ground truth VCF for haplotype 2 (`T_X2_<id>.truth.vcf`). If there is no entry in the ground truth VCF for the haplotype in question the four fields will all contain a '.', indicating that the base at this locus corresponds to that found in the donors personalised reference. False positive somatic variant calls, where they occur, are generally as a result of sequence error, ie., neither haplotype field set in the Ground truth map file record in question is annotated as PASS or MASKED (indicating no variant was spiked-in) and one or both record a sequence error at that locus. It should be noted that Mutect2 false positive calls are rare (particularly in the toy example included, given the small target area in question). Very rarely a false positive call may arise from an alignment artefact. Usually in that case both sets of haplotype fields (4 to 7 and 8 to 11) all contain '.', the variant is passed by Mutect2 (field 3 == PASS)  and the variant in question does not appear in the associated germline VCF. The alignment artefact is confirmed by examining the pileup in the realigned BAM and tracking the QNAME of the reads supporting the "variant" to their true location in the original phased set of BAMs created after spike-in.

![gtMapFormat](https://user-images.githubusercontent.com/63290680/220073276-873c47ec-b6aa-49f8-bbc9-ca19bb8e758f.png)
*Sample output from the Ground truth map file.*

False negatives encountered using Mutect2 significantly outweigh the number of false positives. The ground truth map file allows for easy identification of filtered false negatives in Mutect2's output. False negative entries look very similar to true positive entries. One of the two haplotype field sets will be annotated as PASS or MASKED (indicating a variant was spiked-in). Field 3 however, rather than being set to "PASS" will instead contain a Mutect2 filter reason (for example, "weak_evidence"). The description of each of the Mutect2 filter fields is listed in the metadata near the top of the Mutect2 filtered VCF.

False positive and negative calls can be easily filtered from the ground truth map file using the awk command. For example, to extract false negative calls,
```
$ awk '{if($3 !~ /PASS/) if($6 ~ /PASS/ || $6 ~ /MASK/ || $10 ~ /PASS/ || $10 ~ /MASK/) print $0}' gtMapper.hap.ref

chr19:10145096	0.579	clustered_events	chr19:10144546	.	.	.	chr19:10144496	G>T	PASS	1
chr19:10170621	0.054	normal_artifact;weak_evidence	chr19:10170071	G>A	PASS	0.049809	chr19:10170021	.	.	.
chr19:10228685	0.052	weak_evidence	chr19:10228133	A>T	PASS	0.128522	chr19:10228082	A>G	SEQ_ERROR	0.0454545
```
To extract false positive calls (if there are any),
```
$ awk '{if($3 ~ /PASS/) if($6 !~ /PASS/ && $6 !~ /MASK/ && $10 !~ /PASS/ && $10 !~ /MASK/) print $0}' gtMapper.hap.ref

chr21:15827077	0.165	PASS	chr21:15826798	A>C	SEQ_ERROR	0.285714	chr21:15826328	.	.	.
chr21:37472816	0.231	PASS	chr21:37471383	.	.	.	chr21:37471094	T>G	SEQ_ERROR	0.25
chr21:37511939	0.124	PASS	chr21:37510520	G>C	SEQ_ERROR	0.153846	chr21:37510231	.	.	.
```
Foe example, to extract possible false positives due to mapping artefacts (very rare),
```
$ awk '{ if($3=="PASS") if($4==".") if($5==".") if($6==".")if($8==".") if($9==".") if($10==".") print $0}' gtMapper.hap.ref
```
The summary.txt file, also created by the prepGtMap.bash script and discussed in the next section collates informations similar to this in a brief report output.

### The summary file
The summary is a text file that records a log of the standard output from the prepGtMap.bash command. It provides an overview in plain English on how many true somatic variants were simulated, how many were identified by the caller, how many were filtered incorrectly and how many were missed completely (ie., how many left no record whatsoever in the callers output VCF). It also provides a breakdown of any false positives.

### Filtered false negatives and VAF pdf plots.
This Filtered false negatives pie chart breaks down the contributions of each Mutect2 filter to false negatives listed in its VCF output (i.e., true somatic variants that were incorrectly filtered by the caller as artefacts). The R source to create the plot is also included in the caller directory (i.e., toyExample/MUTECT). Bear in mind that the false negatives listed as filtered by the caller only represent a fraction of the total number of false negatives which Mutect2 failed to call. At about 50x Mutect2 will ignore most of the low frequency burden without leaving any record (ie., filter reason) in its VCF output. The majority of false negatives not listed as filtered by Mutect2 come from the low frequency end of the spectrum. This is evidenced by the order of magnitude difference in scale of the y-axis between the ground truth variant allele frequency plot (ie., the true frequencies present in the sample) and the variant allele frequency plot as estimated from variant frequencies called by Mutect2 below.


![toyExPlotM2GT](https://user-images.githubusercontent.com/63290680/220210089-a16b8e2b-0003-4238-94eb-92147a685150.png)
*Selected output from the toy example showing the contribution of each filter to false negatives and a set of VAF plots. The first VAF plot shows the distribution of the true frequencies present in the sample, the second Mutect2s estimation of those frequencies from the subset of variants it detected.*

### Limitations of the toy example
As the name suggests the example included with this framework is a toy example. It can be used as a base for constructing more comprehensive simulations, however the researcher should be aware of the following limitations and if required address them to suit their application. All of the issues mentioned below can be easily addressed, however a number of them will likely result in an increased run-time.

- The toy example simulations involve a diploid tumour. It is possible to simulate more complex duplication / rearrangement events with this framework however that also makes the simulation more complex. It means you now need to construct a custom tumour reference genome which contains additional contigs to cover duplications and also, perhaps, with some existing contigs removed (if the cancer cell has lost chromosomes etc.). Using a custom tumour reference also changes the format of the ground truth map file. Rather than there being a one-to-two map (one variant locus in the caller output VCF maps to 2 other ground truth loci, one in each haplotype), for duplications your caller VCF locus now maps to 3 or more ground truth loci (1 in each haplotype, and one/more in each duplicated chromosome(s)/region(s)). You obviously need to plan your spike-in config much more carefully. You need to choose the genomic range and haplotype on which a duplication occurs. You need to take into account which somatic variants the duplicated contig picked up before the duplication event and which it picked up after. I hope to give an example of this more complex simulation at some point in the not too distant future and to simplify the framework to make it easier to implement. Please let me know if such a simulation would interest you and in general terms, what your tumour genome of interest would look like.
- The toy example simulates malignant single base substitutions (SBS) only in the tumour BAM. Germline SBS (ie., SNPs) and indels are also simulated in the personalised reference (& therefore in all BAM files produced from it) but malignant indels are not simulated in the tumour BAM. The framework does not currently provide an easy way to do this. It currently requires constructing a custom tumour reference (see previous point). It will be dealt with in a future release however. Similarly the framework does not guarantee phasing of SBS that are located within a read length of each other. They are obviously phased in one haplotype or the other, however both, or just one of variants that close together could end up in the same or separate reads at the pileup locus. This hasn't been an issue for our research to date. It could be addressed with a custom tumour reference but that is not straightforward. It will be addressed in a future release.
- The part of the toy example pipeline that runs Mutect2 (somatic caller) takes some short cuts in the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-). Cross-sample contamination and Base Quality Score Recalibration (BQSR) stages are not run (as our data is simulated without contamination or systematic biases). The functionator is not run (given these are randomly selected target loci, we are not interested in functionator output). The **chr19 genome** included with the toy simulation is a subset of the GRCh38 standard reference. It was small enough to include on the github page. A small genome also speeds up the calling process. In the real world you are obviously not going to be sure all your reads come from chr19 so, if you haven't already, download the standard GRCh38 / latest reference plus indexes and use those in your simulations (i.e., `wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834; tar zxvf 254f697d-310d-4d7d-a27b-27fbf767a834; bwa index -a bwtsw GRCh38.d1.vd1.fa` ). The germline resource used with Mutect2 was also subsetted in the same way (PASS variants from chr19 only) and for the same reason. Feel free to use the full germline resource. These changes shouldn't affect the simulation output in any significant way (other than increasing runtime), however feel free to run all stages in the best practices pipeline with the complete reference and resource files if you wish.
- The **target area** is very small (500KB). It works well for a demo. However, most gene panels are probably around four times that (1 or 2 MB), an exons only exome is around 35MB. You can by all means increase it but as you move past gene panel BEDs and get to exome capture BEDs e.t.c. you will probably need to run this on a cluster.
- The **coverage of 50X** is low for somatic mutation calling. It makes for a quicker simulation but I think, in general, 100X should be the minimum for  calling somatic mutations. Consider increasing it in your simulations.
- The **target BED** used in the toy example is over simplistic (just a contiguous 500KB chunk of chromosome 19). You can of course make your target regions as numerous and complex as you require. Bear in mind however that in these simulations, as in real sequencing data, coverage at the beginning of a capture region does not start at 100X (or whatever) and continue at 100x until the last locus in that target. It ramps up from the start, levels off in the middle and tails off near the end. If your overall target is fragmented into many regions (like for example a merged BED of exons) you are going to end up with less than the average depth of coverage you requested (in `generatePhasedBams.bash`). To work around this, consider adding a flanking region at the start and end of each range (the last optional argument in `createPersonalisedMaskedTarget.sh`). This will pad it out a bit, causing nearby ranges to merge and giving you more even coverage. You may also want to increase the overall coverage a bit to compensate, or use a combination of both approaches. When it comes to dealing with this issue, just like with target capture of real data, the choice is yours.
- The **burden** used in the toy example is very high (TMB = 850 mut/MB). This explains the number of candidate variants filtered with "clustered_events". While there are tumour samples with this level of burden, they are rare. The reason a high burden was chosen was to compensate for the small target area and low depth of coverage of the toy simulation. A lower burden on a larger target and higher depth should still yield interesting results. For example a 100X exome target, 50 base flanking, with 8 mut/MB would probably yield about 5 or 10 times as many filtered false negatives as the toy example. Set the level of burden according to your requirements / model of tumour evolution.
- To keep things simple the **tumour purity** (a.k.a. cellularity) in the toy example is set at 100%. To simulate a different tumour purity simply multiply the spike-in frequencies (in .spike config file) by the tumour purity, expressed as a fraction. For example, for 75% tumour purity, `awk -v tumourPurity=0.75 '{ print $1"\t"$2"\t"$3"\t"($4*tumourPurity) }' chr19_500KB.h1.854mutMB.spike > chr19_500KB.h1.854mutMB.0.75purity.spike # Repeat for h2` Once you have created new configs for both haplotypes, re-run the spike-in, realignment and calling stages of the pipeline to complete the simulation. Please **do not try to adjust by multiplying the allele frequencies estimated by the caller by the purity** as to do so would give you a false and inflated impression of the caller's capability to detect low frequency somatic variants. Instead adjust the spike-in frequencies and re-run the required simulation stages as stated.
- For more extensive simulations using this framework have a look at the github detailing **methods used in the publication** associated with this framework (link to follow..) or contact BrianOSullivan@yahoo.com .


## List of stochasticSim commands and their formats

Some of the more common commands used when setting up a simulation with this framework are shown below. **Please note that all VCF, BED, BAM and spike-in config files must be sorted in chromosome and coordinate order before using them with any of these commands. BAM files must also be indexed.** For a simple example of a script to run a simulation see stochasticSim-main/toyExample/run.bash .

### createPersonalisedMaskedTarget.bash
**SYNOPSIS**

The format is,

```
$ createPersonalisedMaskedTarget.bash <ID prefix> \
                                         <sex, m|f>
                                             <standard reference genome>
                                                 <phased germline VCF>
                                                     <target BED>
                                                         <flanking region, if required>

```
where,
| Argument | Description |
| --- | --- |
| ID prefix | An output file prefix string that will be used as an identifier for all output filenames.|
| sex, m\|f | Sex of the donor. Y chromosome included in simulation if male. Otherwise both phased genome fasta files contain only an X.|
| standard reference genome | Location of the standard reference genome|
| phased germline VCF | Location of donor germline VCF containing SNV and indel information.|
| target BED | Location of BED file specifying target region|
| flanking region | If specified, additional number of bases extending on either side of each range in the target BED file to be included in the simulation. May be useful to simulate off-target capture of just to increase the target region size if required.|

**DESCRIPTION**

Create a phased, personalised reference genome that is modified to reflect target region of interest. Usually executed as the first stage in the simulation process. All locations outside the target area will have their contents replaced with 'N' bases to avoid reads being simulated from these regions. In addition to the pair of phased fasta reference files this command also output a number of other files for use in subsequent steps.   
    
For example,
```
../bin/createPersonalisedMaskedTarget.bash HG00110 F GRCh38.d1.vd1.chr21.fa HG00110.chr21.vcf chr21.exons_ranges.bed
```
This command will output the following files,

| Filename | Description |
| --- | --- |
| X1_\<ID prefix\>.fa | reference assembly for haplotype 1|
| X\|Y2_\<ID prefix\>.fa |  reference assembly for haplotype 2|
| X1_\<ID prefix>\.\<target BED\>.fa | reference assembly for haplotype 1, masked to only include target region|
| X\|Y2_\<ID prefix\>.\<target BED\>.fa | reference assembly for haplotype 2, masked to only include target region|
| X1_\<ID prefix\>.\<target BED\>.merged.bed.gz | Target region lifted over to haplotype 1 reference assembly|
| X\|Y2_\<ID prefix\>.\<target_BED\>.merged.bed.gz | Target region lifted over to haplotype 2 reference assembly|
| liftover_X1_\<ID prefix\>.txt | A liftover file, enabling mapping between genomic locations in the standard reference and those in reference assembly for haplotype 1|
| liftover_X1_\<ID prefix\>.txt | A liftover file, enabling mapping between genomic locations in the standard reference and those in reference assembly for haplotype 2|

### condenseLift
**SYNOPSIS**
This (binary) converts the liftover*.txt files created by the previous command to a much smaller and compact format for use with subsequent commands. Once the conversion is complete the original liftover*.txt may be deleted. Usually executed as the second step in the simulation, directly after `createPersonalisedMaskedTarget.bash`.

The format is,

```
$ ../bin/condenseLift \<large format liftover file\>

```
where,
| Argument | Description |
| --- | --- |
| liftover_X\|Y1\|2_HG00110.txt | large format liftover file |

**DESCRIPTION**
  
    
For example,
```
$ ../bin/condenseLift liftover_X1_HG00110.txt > liftover_X1_HG00110.condense.txt
```
This command outputs to the standard out which may be redirected to a file as shown in the example above.


### generatePhasedBams.bash
**SYNOPSIS**

The format is,

```
$ generatePhasedBams.bash <haplotype 1 reference>
                              <haplotype 2 reference>
                                  <read length>
                                     <depth>
                                         <mean fragment length>
                                             <user suffix>

```
where,
| Argument | Description |
| --- | --- |
| haplotype 1 reference | The personalised target genome, as created in the previous step. |
| haplotype 2 reference | The second (phased) personalised target genome, as created in the previous step. |
| read length | Length of reads to be simulated. |
| Coverage | Coverage argument passed within parenthesis to ART. May be fold ("-f \<fold coverage\>") or read cound ("-c \<number of reads\>"). |
| mean fragment length | Average insert length of data to be simulated. |
| user suffix | An output file prefix string that will be used as an identifier for all output filenames. |

**DESCRIPTION**

This simulates phased sets of tumour and normal BAMs plus index files that will be used as a base to spike in the required somatic distribution. Usually executed as the third step in the simulation process.
    
For example,
```
$ ../bin/generatePhasedBams.bash X1_HG00110.chr21.exons_ranges.fa X2_HG00110.chr21.exons_ranges.fa 76 "-f 25" 180 "25x_76bp"
```
This command will output the following files,

| Filename | Description |
| --- | --- |
| T_X1_\<ID prefix\>.bam | Haplotype 1 pre-Tumour BAM file |
| T_X1_\<ID prefix\>.bam.bai | Haplotype 1 pre-Tumour BAM index |
| T_X2_\<ID prefix\>.bam | Haplotype 1 pre-Tumour BAM file |
| T_X2_\<ID prefix\>.bam.bai | Haplotype 1 pre-Tumour BAM index |
| N_X1_\<ID prefix\>.bam | Haplotype 1 Normal BAM file |
| N_X1_\<ID prefix\>.bam.bai | Haplotype 1 Normal BAM index |
| N_X2_\<ID prefix\>.bam | Haplotype 1 Normal BAM file |
| N_X2_\<ID prefix\>.bam.bai | Haplotype 1 Normal BAM index |

### spikeIn.bash
**SYNOPSIS**

The format is,

```
$ spikeIn.bash  <haplotype 1 pre-tumour BAM> \
                  <haplotype 2 pre-tumour BAM>  \
                     <haplotype 1 reference> \
                        <haplotype 2 reference> \
                           <haplotype 1 tumour somatics list> \
                              <haplotype 2 tumour somatics list>
```
where,
| Argument | Description |
| --- | --- |
| haplotype 1 pre-tumour BAM | A pre-tumour BAM simulated from haplotype 1 reference. |
| haplotype 1 pre-tumour BAM | A pre-tumour BAM simulated from haplotype 2 reference. |
| haplotype 1 reference | Fasta file detailing the first set of chromosomes in the personalised, diploid reference assembly. |
| haplotype 2 reference | Fasta file detailing the second set. |
| haplotype 1 tumour somatics list | A tab seperated file listing the genomic location and types of SNVs to be simulated in hap. 1. |
| haplotype 2 tumour somatics list | Equivalent tab seperated file for hap. 2. |

**DESCRIPTION**

This command spikes in the required burden to pre-tumour BAM files. In real genomic sequencing data the observed number of reads containing the alternate allele is a random variable. This is also the case when spiking in variants with the simulation framework. With spikeIn.bash (& the binary stochasticSpike it calls), the probability of any given read at the target pileup picking up the synthetic alternate allele is taken as the true alternate allele frequency. By accounting for the random nature of genomic sequencing, our simulations more accurately reflect the frequency profile of real tumour genomic sequencing data. This command is usually executed as the fourth step in the simulation process after the pre-tumour BAMs have been created.

For example,
```
$ ../bin/spikeIn.bash T_X1_HG00110.chr21.exons_ranges_25x_76bp.bam \
                        T_X2_HG00110.chr21.exons_ranges_25x_76bp.bam \
                           X1_HG00110.chr21.exons_ranges.fa X2_HG00110.chr21.exons_ranges.fa \
                              chr21.h1.TMB_743mutMb.NE.spike \
                                 chr21.h2.TMB_743mutMb.NE.spike
```
Note that input BAM files must be indexed. This command will output the following files,

| Filename | Description |
| --- | --- |
| haplotype 1 tumour BAM | A tumour BAM (from haplotype 1) containing the required burden. |
| haplotype 2 tumour BAM | A tumour BAM (from haplotype 2) containing the required burden. |
| tumour hap 1 index | Tumour BAM (from haplotype 1) index. |
| tumour hap 2 index| Tumour BAM (from haplotype 2) index. |

### realignAndMerge.bash
**SYNOPSIS**

The format is,

```
$ realignAndMerge.bash  <haplotype 1 BAM> \
                         <haplotype 2 BAM>  \
                           <reference> \
                             <output BAM prefix>

```
where,
| Argument | Description |
| --- | --- |
| haplotype 1 BAM | BAM file with hap. 1 burden, aligned to personalised reference for hap. 1 |
| haplotype 2 BAM | BAM file with hap. 2 burden, aligned to personalised reference for hap. 2 |
| reference | the standard reference (i.e., GRCh38.d1.vd1.fa) |
| output BAM prefix | What you want to call the realigned, combined output BAM |

**DESCRIPTION**

This command is used to realign and merge a tumour or normal phased BAM pair. It is usually the final stage in the simulation process before somatic variant calling.
    
For example,
```
$ ../bin/realignAndMerge.bash T_X1_HG00110.exons_ranges_25x_76bp.spike.bam  \
                                T_X2_HG00110.exons_ranges_25x_76bp.spike.bam \
                                  GRCh38.d1.vd1.fa \
                                    T_HG00110.exons_ranges_50x_76bp
```
This command will output the following file,

| Filename | Description |
| --- | --- |
| \<output BAM prefix\>.bam | Tumour BAM with burden spiked in and aligned against standard reference |
| \<output BAM prefix\>.bam.bai | Corresponding tumour BAM index |

### xxx.bash
**SYNOPSIS**

The format is,

```
$ 

```
where,
| Argument | Description |
| --- | --- |

**DESCRIPTION**
  
    
For example,
```
$ ../bin/
```
This command will output the following files,

| Filename | Description |
| --- | --- |
