# Stochastic simulation framework<!-- omit in toc -->

*A set of software tools providing comprehensive and realistic simulation of tumour genomic sequencing data.*

Contact BrianOSullivan@yahoo.com with questions.

## Table of contents<!-- omit in toc -->
- [System requirements]
- [Installation]
- [Getting started]

## System requirements

This simulation framework runs on a 64bit, Mac OS X or GNU/Linux platform.
If this does not cover your system this script will not work.
If you need a version for your system let me know (BrianOSullivan@yahoo.com).
If you want to run Mutect2 as your somatic caller (as in toyExample included)
you will also need Java 8 / JRE or SDK 1.8 (check your version with "java -version").
The script "prepGtMap.bash" script included plots the contribution of each GATK filter
to the number of false negatives in the output VCF post variant calling.
For this script to run you must have R (script, /usr/bin/Rscript)
installed on your system together with libraries ggplot2, ggrepel, vctrs and forcats.

The total disc space taken up by the framework, once all associated tools
(samtools, bwa etc.) have been locally installed and the toy example simulation run
is about 1.7G. Aprox. 0.7G of this is taken up by BAM files created during simulation.
These files may be deleted once the somatic variant caller has run
(unless you need to look into alignment artefacts in caller output).

This install builds the following publically available tools for use with this framework

      * bwa 0.7.17.
      * samtools 1.13
      * bedtools 2.29.2

You may if you wish use other versions of these tools you already have installed
on your system. They will probably work, however as I have not tested it, I can not
guarantee it. Links to all scripts and binaries used in this simulation framework
are included under the top level bin directory (stochasticSimFramework/bin).
You can copy them into another directory in your path and run them from there
if you prefer.

This script will also download ART read simulator (mountrainier2016) binaries
and the gatk-package-4.2.2.0 jar file (containing MUTECT2 somatic variant caller).


## Installation

Before running this script, run the following commands...
Make a directory into which the simulation framework and the set of tools on which it depends will be installed.
```
mkdir stochasticSimFramework
cd stochasticSimFramework
```

Now download the simulation framework, either latest release of clone head of tree from gitHub.
```
# tar zxvf ../../stochasticSim.tgz
```

Next cd to the stochasticSim directory. Look at the installation script to find out what it will install.
Finally run the install script.
```
cd stochasticSim
more ./install.bash
./install.bash
```

## Getting started

