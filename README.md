# Stochastic simulation framework<!-- omit in toc -->

*A set of software tools providing comprehensive and realistic simulation of tumour genomic sequencing data.*

Contact BrianOSullivan@yahoo.com with questions.

## Table of contents<!-- omit in toc -->
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Getting started](#getting-started)

## System requirements

This simulation framework runs on a 64bit, Mac OS X or GNU/Linux. If this does not cover your system it will not work. If you need a version for a different platform let me know (BrianOSullivan@yahoo.com). If you want to run Mutect2 as your somatic caller (as in toyExample included) you will also need Java 8 / JRE or SDK 1.8 installed (check your version with "java -version"). The "prepGtMap.bash'' script included which plots the contribution of each GATK filter to the number of false negatives in the VCF output also requires R (script, /usr/bin/Rscript) together with libraries ggplot2, ggrepel, vctrs and forcats (see https://cran.r-project.org if you do not already have R installed).

The total disc space taken up by the framework, once all associated tools (samtools, bwa etc.) are installed locally and the toy example simulation run is about 1.7G. Approximately 0.7G of this is taken up by BAM files created during simulation. These files may be deleted once the somatic variant caller has run (unless you need to look into alignment artefacts/issues in caller output, in that case you will need them).

This install also builds the following publically available tools for use with this framework.

    bwa >0.7.17.
    samtools 1.13
    bedtools 2.29.2

They are installed locally, in the same directory in which stochasticSim was cloned and neither the install or simulation process require root access to run. The installation process will not modify your path or cause any conflict with any previous installation of these tools you may have on your system. You may if you wish use other versions of these tools with stochasticSim if you have them installed. They will probably work, however as I have not tested it, I can not guarantee it. Links to all scripts and binaries used in this simulation framework are included under the top level bin directory (stochasticSimFramework/bin). You can copy them into another directory in your path and run them from there if you prefer. The install script uses curl to download these dependencies. You must have curl installed on your system or the install will not work (ie., for ubuntu linux, apt install curl). There are issues when building bwa aligner version 0.7.17 with [some gcc versions (> 10)](https://github.com/lh3/bwa/pull/385) so we use a clone of the bwa main branch in the install which contains a fix.

To build stochasticSim and associated tool dependencies you will need a basic set of developer resources set up such as GNU make, C / C++ compiler, autoconf e.t.c. On many installations these will already be in place (particularly if they have been used to build tools previously). If not however you will need to install them (with for example apt-get install on ubuntu linux, or homebrew / the appropriate developer install for Mac). These installations will require root access. Samtools (& HTSlib) depend on a number of libraries. Again, it is likely these are already on you system. On ubuntu linux you can check with,



Finally, the install script will also download ART read simulator (mountrainier2016) binaries and the gatk-package-4.2.2.0 jar file (containing MUTECT2 somatic variant caller, requires Java 8 / JRE or SDK 1.8).



## Installation

Before running this script, run the following commands.
Make a directory into which the simulation framework and the set of tools on which it depends will be installed. Enter that directory
```
cd <your download target directory>
mkdir stochasticSimFramework
cd stochasticSimFramework
```

Now download and unpack the simulation framework, either latest release or main branch from gitHub, for example using main,
```
unzip <download path>/stochasticSim-main.zip
```

Enter the stochasticSim directory. Look at the installation script to find out what it will install.
Run the install script, checking for any errors in the output and if required rectify any issues that may have occurred during installation.
```
cd stochasticSim
more ./install.bash
./install.bash
```

## Getting started
Once the install has completed successfully you can test it out on the toy example included.
First source the 'tool.path' file to ensure the simulation uses
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
Look at <install path>/stochasticSimFramework/stochasticSim-main/toyExample/MUTECT for the simulations results.
Along with the variant caller output, this directory will contain the following files

| Filename | Description |
| --- | --- |
| gtMapper.hap.ref | A tab seperated table mapping all entries in the caller filtered VCF output to their ground truth values in each haplotype. |
| plotGtPieChart.R/pdf | Source and associated pdf output to plot a pie chart of all caller filtered false negatives. |
| summary.txt | An overall breakdown of where caller false positives/negatives occurred in this simulation and why. |

