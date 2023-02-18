# Stochastic simulation framework<!-- omit in toc -->

*A set of software tools providing comprehensive and realistic simulation of tumour genomic sequencing data.*

Contact BrianOSullivan@yahoo.com with questions.

## Table of contents<!-- omit in toc -->
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Getting started](#getting-started)

## System requirements

This simulation framework runs on a 64bit, Mac OS X or GNU/Linux. If this does not cover your system it will not work. If you need a version for a different platform let me know (BrianOSullivan@yahoo.com). If you want to run Mutect2 as your somatic caller (as in toyExample included) you will also need Java 8 / JRE or SDK 1.8 installed (check your version with `java -version`). The `prepGtMap.bash` script included which plots the contribution of each GATK filter to the number of false negatives in the VCF output also requires `R` (script, `/usr/bin/Rscript`) together with libraries `ggplot2`, `ggrepel`, `vctrs` and `forcats` (see https://cran.r-project.org if you do not already have R installed).

The total disc space taken up by the framework, once all associated tools (samtools, bwa etc.) are installed locally and the toy example simulation run is about 1.7G. Approximately 0.7G of this is taken up by BAM files created during simulation. These files may be deleted once the somatic variant caller has run (unless you need to look into alignment artefacts/issues in caller output, in that case you will need them).

This install also builds the following publically available tools for use with this framework.

    bwa >0.7.17.
    samtools 1.13
    bedtools 2.29.2

They are installed locally, in the same directory in which `stochasticSim` was cloned and neither the install or simulation process require root access to run. The installation process will not modify your path or cause any conflict with any previous installation of these tools you may have on your system. You may if you wish use other versions of these tools with `stochasticSim` if you have them installed. They will probably work, however as I have not tested it, I can not guarantee it. Links to all scripts and binaries used in this simulation framework are included under the top level bin directory (`stochasticSimFramework/bin`). You can copy them into another directory in your path and run them from there if you prefer. The install script uses curl to download these dependencies. You must have curl installed on your system or the install will not work (ie., for ubuntu linux, `apt install curl`). There are issues when building bwa aligner version 0.7.17 with [some gcc versions (> 10)](https://github.com/lh3/bwa/pull/385) so we use a clone of the bwa main branch in the install which contains a fix.

To build `stochasticSim` and associated tool dependencies you will need a basic set of developer resources set up such as GNU make, C / C++ compiler, autoconf e.t.c. On many installations these will already be in place (particularly if they have been used to build tools previously). If not however you will need to install them (with for example `apt-get install` on ubuntu linux, or homebrew / the appropriate developer install for Mac). These installations will require root access. Samtools (& HTSlib) depend on a number of libraries. Again, it is likely these are already on you system. On ubuntu linux for example you can check with,

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

Two ground truth VCFs files (one for each haplotype, *.truth.vcf) are created during spike-in against perfectly aligned, phased data. The filter field in each record indicates the true reason behind all non-reference loci in the data.

| Filter field | Description |
| --- | --- |
| SEQ_ERROR | The locus does not contain an alternative allele. The non reference base(s) is a result of sequence error during simulation |
| PASS | The locus contains a somatic variant and all reads supporting that variant were successfully spiked into the data |
| MASKED | The locus contains a somatic variant however one or more supporting reads in the data were masked by sequence error |
| NO_COVERAGE | The locus contains a somatic variant however no coverage was obtained at this locus during simulation (DP=0)|
| UNDETECTED | The locus contains a somatic variant however no reads containing the alternative allele were detected during the simulation (AD for the spiked-in allele=0) |

The ground truth VCFs are created prior to somatic variant calling, when the simulated data is perfectly aligned against the phased personalised donor reference. Issues that arise after the data is realigned against the standard reference and used as input to a somatic variant caller are recorded in the ground truth map file (`gtMapper.hap.ref`).

### Ground truth map file

This file maps each record in the (filtered VCF) somatic variant caller output to its corresponding ground truth in each haplotype (recorded in the ground truth VCFs). The first 2 columns of each record in the ground truth map show the genomic coordinates (ie., chromosome and locus, 1 base) of the variant as it appears in the caller output VCF. These coordinates correspond to its aligned location in the standard reference genome.

