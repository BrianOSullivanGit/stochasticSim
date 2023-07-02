#!/bin/bash

# Install script for stochastic simulation framework.

# Requirements: 64bit, Mac OS X or GNU/Linux platform
#               If this does not cover your system this script will not work.
#               If you need a version for your system let me know (BrianOSullivan@yahoo.com)
#               If you want to run Mutect2 as your somatic caller (as in toyExample included)
#               you will also need Java 8 / JRE or SDK 1.8 (check your version with "java -version").
#               The script "prepGtMap.bash" which is run after (Mutect2) variant calling
#               creates an R script. This script plots the contribution of each GATK filter
#               to the number of false negatives in the output VCF.
#               For this script to run you must have R (script, /usr/bin/Rscript)
#               installed on your system together with libraries ggplot2, ggrepel, vctrs and forcats.
#
#               The total disc space taken up by the framework, once all associated tools
#               (samtools, bwa etc.) have been locally installed and the toy example simulation run
#               is about 1.7G. Aprox. 0.7G of this is taken up by BAM files created during simulation.
#               These files may be deleted once the somatic variant caller has run
#               (unless you need to look into alignment artefacts in caller output).
#
#               This install builds the following publically available tools for use in this framework
#
#                      bwa 0.7.17.
#                      samtools 1.13
#                      bedtools 2.29.2
#
#               You may if you wish use other versions of these tools you already have installed
#               on your system. They will probably work, however as I have not tested it, I can not
#               guarantee it. Links to all scripts and binaries used in this simulation framework
#               are included under the top level bin directory (stochasticSimFramework/bin).
#               You can copy them into another directory in your path and run them from there
#               if you prefer.
#               
#               This script will also download ART read simulator (mountrainier2016) binaries
#               and the gatk-package-4.2.2.0 jar file (containing MUTECT2 somatic variant caller).	
	


# Before running this script, run the following commands...

# Make a directory into which the simulation framework and the set
# of tools on which it depends will be installed.

# mkdir stochasticSimFramework
# cd stochasticSimFramework

# Now download the simulation framework, either latest release of clone head of tree from gitHub.
# tar zxvf ../../stochasticSim.tgz

# Finally cd to the stochasticSim directory and run the install script.
# cd stochasticSim
# ./install.bash

# Install starts here......

# We assume you've followed instructions and cd'd to the stochastimSim directory to run this script.
stochasticSimDirName=`pwd|sed 's/.*\///1'`

# Install other tool dependancies in parent directory.
cd ../
mkdir bin

date | tr '\012' ':'
echo " Setting up framework dependancies."

# Download & build required tools.

# Bwa
date | tr '\012' ':'
echo " Build bwa head."
git clone https://github.com/lh3/bwa.git
make -C bwa || { echo -e "\n\033[7mBWA Build failed. Resolve issues before proceeding.\033[0m";exit 1; }
# Leave a link in bin
ln -s ../bwa/bwa bin/bwa
echo "export ALIGNER="`pwd`"/bwa/bwa" >> bin/tool.path


# Samtools
date | tr '\012' ':'
echo " Build samtools 1.13."
curl -L https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 | tar jxvf - || { echo -e "\n\033[7mSAMTOOLS download failed. Resolve issues before proceeding.\033[0m";exit 1; }
cd samtools-1.13
mkdir local_copy_output

# Workaround to install on a system without curses.
if [ -d /usr/include/ncurses ]; then
  # Assume ncurses is there, install as normal
  ./configure --prefix=${PWD}/local_copy_output
else
  # User does not want or is unable to install ncurses, try without,
  ./configure --prefix=${PWD}/local_copy_output --without-curses
fi

make || { echo -e "\n\033[7mSAMTOOLS Build failed. Resolve issues before proceeding.\033[0m";exit 1; }
make install || { echo -e "\n\033[7mSAMTOOLS install failed. Resolve issues before proceeding.\033[0m";exit 1; }
# This is needed by sim. framework that links htslib.
SAMTOOLS_BUILD_PATH=${PWD}
HTSLIB_VERSION=`ls -d htslib-* | sed 's/.*-//'`
export SAMTOOLS_BUILD_PATH
export HTSLIB_VERSION
cd ../
ln -s ../samtools-1.13/local_copy_output/bin/samtools bin/samtools
echo "export SAMTOOLS="`pwd`"/samtools-1.13/local_copy_output/bin/samtools" >> bin/tool.path

# ART read simulator
date | tr '\012' ':'
echo " Install ART read simulator (mountrainier2016) binaries."	
# Just install binaries..
if [ "$(uname)" == "Darwin" ]; then
    # Install under Mac OS X platform
    curl -L https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05macos64.tgz | tar zxvf - || { echo -e "\n\033[7mART download failed. Resolve issues before proceeding.\033[0m";exit 1; } 
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Install under GNU/Linux platform
    curl -L https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz | tar zxvf - || { echo -e "\n\033[7mART download failed. Resolve issues before proceeding.\033[0m";exit 1; } 
else
    echo "Platform \"$(uname)\" not supported...exiting..."
    exit 5
fi

ln -s ../art_bin_MountRainier/art_illumina bin/art_illumina
echo "export READ_SIMULATOR="`pwd`"/art_bin_MountRainier/art_illumina" >> bin/tool.path

# Bedtools
date | tr '\012' ':'
echo " Build bedtools 2.29.2"	
curl -L https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz | tar zxvf - || { echo -e "\n\033[7mBEDTOOLS download failed. Resolve issues before proceeding.\033[0m";exit 1; } 
cd bedtools2/
# Fix for bedtools install issue with some verion of Mac (clang) compiler.
# There may also be a problem with python "not found" (ie., a 2.7 -> 3.0 thing)
# Consider installing python-is-python3 if you're on ubuntu.
# 
# Putting in a temp fix for now.

if [ "$(uname)" == "Darwin" ]; then
    # Install under Mac OS X platform
    # Workaround for compiler problem
    mv ./src/utils/gzstream/version ./src/utils/gzstream/version.txt
fi

python --version 2>&1 > /dev/null || { echo -e "\n\033[7mWARNING Can't find a version of python available to this script (is python3 linked to python?). This may cause problems.\033[0m"; }

make || { echo -e "\n\033[7mWARNING BEDTOOLS Build failed. Trying to proceed with bedtools binary that was created. Please check and resolve any issues.\033[0m"; }

if [ -x bin/bedtools ]; then
  # OK it built the bedtools binary at least, try to work with that.
  cd ../
  ln -s ../bedtools2/bin/bedtools bin/bedtools
  echo "export BEDTOOLS="`pwd`"/bedtools2/bin/bedtools" >> bin/tool.path
else
  echo -e "\n\033[7mFailed to build BEDTOOLS binary. Please resolve issues before proceeding.\033[0m"
  exit 1
fi

# Datamash
curl -L http://ftp.gnu.org/gnu/datamash/datamash-1.3.tar.gz | tar zxvf - || { echo -e "\n\033[7mDATAMASH download failed. Resolve issues before proceeding.\033[0m";exit 1; } 
cd datamash-1.3
./configure
make || { echo -e "\n\033[7mDATAMASH Build failed. Resolve issues before proceeding.\033[0m";exit 1; }
make check || { echo -e "\n\033[7mWARNING DATAMASH Build check failed. Trying to proceed with datamash binary that was created. Please check and resolve any issues.\033[0m"; }

if [ -x datamash ]; then
  # OK it built the datamash binary at least, try to work with that.
  cd ../
  ln -s ../datamash-1.3/datamash bin/datamash
  echo "export DATAMASH="`pwd`"/datamash-1.3/datamash" >> bin/tool.path
else
  echo -e "\n\033[7mFailed to build DATAMASH binary. Please resolve issues before proceeding.\033[0m"
  exit 1
fi


# Now build stochastic sim. framework
date | tr '\012' ':'
echo " Build ${stochasticSimDirName}"	
cd ${stochasticSimDirName}
make all || { echo -e "\n\033[7mSTOCHASTICSIM Build failed. Resolve issues before proceeding.\033[0m";exit 1; }
cd ../
ln -s ../${stochasticSimDirName}/bin/createDonorGenome bin/createDonorGenome
ln -s ../${stochasticSimDirName}/bin/liftover bin/liftover
ln -s ../${stochasticSimDirName}/bin/stochasticSpike bin/stochasticSpike
ln -s ../${stochasticSimDirName}/bin/targetRef bin/targetRef
ln -s ../${stochasticSimDirName}/bin/createPersonalisedMaskedTarget.bash bin/createPersonalisedMaskedTarget.bash
ln -s ../${stochasticSimDirName}/bin/generatePhasedBams.bash bin/generatePhasedBams.bash
ln -s ../${stochasticSimDirName}/bin/prepGtMap.bash bin/prepGtMap.bash
ln -s ../${stochasticSimDirName}/bin/realignAndMerge.bash bin/realignAndMerge.bash
ln -s ../${stochasticSimDirName}/bin/spikeIn.bash bin/spikeIn.bash
ln -s ../${stochasticSimDirName}/bin/vcfAntex bin/vcfAntex
ln -s ../${stochasticSimDirName}/bin/tncSpike bin/tncSpike
ln -s ../${stochasticSimDirName}/bin/tncCountsProfile bin/tncCountsProfile

echo "export CREATEDONORGENOME="`pwd`"/${stochasticSimDirName}/bin/createDonorGenome" >> bin/tool.path
echo "export LIFTOVER="`pwd`"/${stochasticSimDirName}/bin/liftover" >> bin/tool.path
echo "export STOCHASTIC_SPIKE="`pwd`"/${stochasticSimDirName}/bin/stochasticSpike" >> bin/tool.path
echo "export TARGETREF="`pwd`"/${stochasticSimDirName}/bin/targetRef" >> bin/tool.path
echo "export CONDENSE_LIFT="`pwd`"/${stochasticSimDirName}/bin/condenseLift" >> bin/tool.path
echo "export GTMAPPER="`pwd`"/${stochasticSimDirName}/bin/gtMapper" >> bin/tool.path
echo "export TWOWAYLIFTOVER="`pwd`"/${stochasticSimDirName}/bin/2wayLiftover" >> bin/tool.path
echo "export CREATEPERSONALISEDMASKEDTARGET_BASH="`pwd`"/${stochasticSimDirName}/bin/createPersonalisedMaskedTarget.bash" >> bin/tool.path
echo "export GENERATEPHASEDBAMS_BASH="`pwd`"/${stochasticSimDirName}/bin/generatePhasedBams.bash" >> bin/tool.path
echo "export PREPGTMAP_BASH="`pwd`"/${stochasticSimDirName}/bin/prepGtMap.bash" >> bin/tool.path
echo "export REALIGNANDMERGE_BASH="`pwd`"/${stochasticSimDirName}/bin/realignAndMerge.bash" >> bin/tool.path
echo "export SPIKEIN_BASH="`pwd`"/${stochasticSimDirName}/bin/spikeIn.bash" >> bin/tool.path
echo "export VCFANTEX="`pwd`"/${stochasticSimDirName}/bin/vcfAntex" >> bin/tool.path
echo "export TNCSPIKE="`pwd`"/${stochasticSimDirName}/bin/tncSpike" >> bin/tool.path
echo "export TNCCOUNTSPROFILE="`pwd`"/${stochasticSimDirName}/bin/tncCountsProfile" >> bin/tool.path


# Finally setup the reference, target bed and somatic variant caller that we will use in our toyExample.
date | tr '\012' ':'
echo " Setup ${stochasticSimDirName} toy example."	
cd ${stochasticSimDirName}/toyExample
# Unpack the reference and bed file etc. for this toy example.
tar zxvf GRCh38.d1.vd1.HG00110.chr19.tgz
# Index it.
../../bwa/bwa index -a bwtsw GRCh38.d1.vd1.chr19.fa || { echo -e "\n\033[7mFailed to index GRCh38 chr19 genome. Resolve issues before proceeding.\033[0m";exit 1; }
 
# Setup caller (Mutect2)
date | tr '\012' ':'
echo " Download GATK 4.2.2.0 jar."
curl -L https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip > gatk-4.2.2.0.zip || { echo -e "\n\033[7mGATK JAR download failed. Resolve issues before proceeding.\033[0m";exit 1; } 
unzip gatk-4.2.2.0.zip
rm gatk-4.2.2.0.zip
cd ../../

cp bin/tool.path ${stochasticSimDirName}/bin

# Make sure we have a full  set of tools (there should be X of them TODO!!!! update this!!)

# Check if the links to any of these tools are broken.
ls bin | while read line; do if [ ! -e bin/${line} ]; then echo ${line}; fi; done > $$.prob.lst
if [ -s $$.prob.lst ]; then echo -e "\n\033[7mINSTALL ERROR:\033[0m There was a problem building the following tools,\n"; cat $$.prob.lst;echo -e "\nCheck back through the install output to resolve the issue before proceeding."; rm $$.prob.lst; exit 1; fi
rm $$.prob.lst

#The demo simulation included is based on the allele frequency distribution
#mutant allele frequency spectrum of a diploid tumour expected under a neutral
#evolutionary subclonal model with a clonal point-mass

echo "The stochastic simulation framework setup is complete."
echo "Check the output above and rectify any errors that may have occurred during installation."
echo
echo "Once the install has completed successfully you can test it out on the toy example included."
echo "First source the 'tool.path' file to ensure the simulation uses"
echo "the required versions of associated software tools."
echo
echo " # source ${PWD}/${stochasticSimDirName}/bin/tool.path"
echo
echo "Now, enter the toy example directory and run the simulation,"
echo 
echo " # cd ${PWD}/${stochasticSimDirName}/toyExample"
echo " # ./run.bash 50 chr19_500KB.bed"
echo 
echo "Depending on your hardware this will take about 8 minutes to run."
echo "Look at ${PWD}/${stochasticSimDirName}/toyExample/MUTECT for the simulations results."
echo "Along with the variant caller output, this directory will contain the following files"
echo
echo -e "\033[1mgtMapper.hap.ref:\033[0m     A tab seperated table mapping all entries in the caller filtered VCF output to their\n                      ground truth values in each haplotype."
echo -e "\033[1mplotGtPieChart.R/pdf:\033[0m Source and associated pdf output to plot a pie chart of all caller filtered false negatives."
echo -e "\033[1msummary.txt:\033[0m          An overall breakdown of where caller false positives/negatives occurred in this simulation\n                      and why."
echo
