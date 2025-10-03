#!/bin/bash

#####################################################################
# Copyright (c) 2020 by University of Southern California,
# Institute of Translational Genomics; Author: Yuxin Jin
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#####################################################################

time=`date +%d-%m-%Y-%H-%M`

echo
echo
echo "########################################################################"
echo "########################################################################"
echo "### Starting $0 at $time "
echo "########################################################################"
echo "########################################################################"

# Define variable
runDir=`pwd`
myhostname=`hostname`
sampleConfig=$runDir/SampleID.csv
projList=`grep ^K $sampleConfig | cut -d_ -f1 | sort | uniq | tr "\n" " "`
#projList=`grep ^C $sampleConfig | cut -d_ -f1 | sort | uniq | tr "\n" " "`
slurmHome=/home1/rdagnew/deliveryScripts
baseDir=/project/davidwcr_264/Projects/KGP
runID=`pwd | rev | cut -d/ -f1 | rev`

#Yuxin made some changeds on 08/18/2022. In the future, pipeline release folder will be project based, then run based.

echo "### --Running on $myhostname-- ###"
echo "### Start creating postrun delivery structure on $runDir"


if [[ ! -e ${runDir}/multiqc ]]; then
        echo "### multiqc folder doesn't exist ###"
        echo "### Start create multiqc folder ###"
        mkdir -vp ${runDir}/multiqc
else
	echo "### Found multiqc folder, continue###"
fi


for projID in $projList
do
	# Clarify recipe per project
	echo "########################################################################"
	echo "### Start creating delivery structure for project $projID ###"

	if [ -e echoOut/$projID ] && [ ! -e islipOut/$projID ]
	then
		recipe="echo"
		echo "### Found echoOut, recipe for project $projID should be $recipe"
	elif [ -e islipOut/$projID ] && [ ! -e echoOut/$projID ]
	then
		recipe="islip"
		echo "### Found islipOut, recipe for project $projID should be $recipe"
	elif [ -e islipOut/$projID ] && [ -e echoOut/$projID ]
	then
		recipe="islip echo"
		echo "### Found both islipOut and echoOut, recipe for projects $projID should be $recipe"
	else
		echo "### Something goes wrong!! Please check whether pipeline outdir exists!!! ###"
		echo "### continue !!! ###"
		continue
	fi

	# Create sampleID and base folder for release
	for recipeID in $recipe
	do
		echo "### Creating release folder for ${recipeID}-${projID}"
		mkdir -vp ${runDir}/release/${recipeID}-${projID}
		mkdir -vp ${runDir}/release/${recipeID}-${projID}_RC1
		configPath=${runDir}/release/${recipeID}-${projID}
		ID=`grep "^${projID}" ${sampleConfig} | cut -d_ -f1,2 | sort | uniq | sed "s/${projID}/ID=${projID}/g"`
		kitType=`grep "^${projID}" ${sampleConfig} | cut -d_ -f7 | uniq`
		#demulPath=${baseDir}/${projID}/${runID}/FASTQs
                demulPath=${runDir}/FASTQs

		# Define kit
		case ${kitType} in
			Genome|TSWGL|KHWGSH|KHWGSM )
			kit="Genome"
			;;
			exome|ASS8A|AV6UOS|KHSSEM|A1S4X|ASXE6C|ASXE6S|ASXE6U|KHS6C|KHST2|ASKIL6|KHXG1|K1STX|SCEZEX|KHSSH5|ITPANC|SSXTCP|ASX8NC )
			kit="Exome"
			;;
			longRNA|RNA|SM2V2|TSMRS|KHMRS|OHFSK|KBMRS|OFKHY|TSMRU|TSRAP|A1MRS|KRSOD|TSRGD|NEBKAP|SWI2SU|BGIRNA|UNRNAF|UNRNAS|UNRNAN|BGIRNA|NEBNU2|ILEXPAN )
			kit="RNA"
			;;
			* )
			echo "### Kit type ($kitType) not found.  Please check samplename."
			continue
		esac

		echo "### kitType is $kitType"

		# Create sampleID.txt
		# Check $configPath/${recipeID}-${projID}_sampleID.txt existance

		echo "### Start creating $configPath/${recipeID}-${projID}_sampleID.txt"
		
		if [ -e $configPath/${recipeID}-${projID}_sampleID.txt ]
		then
			echo "### $configPath/${recipeID}-${projID}_sampleID.txt exist"
			echo "### Will skip this step"
		else
			echo "recipe=$recipeID" > $configPath/${recipeID}-${projID}_sampleID.txt
			echo "kitType=$kitType" >> $configPath/${recipeID}-${projID}_sampleID.txt
			echo "dataType=$kit" >> $configPath/${recipeID}-${projID}_sampleID.txt
			echo "project=$projID" >> $configPath/${recipeID}-${projID}_sampleID.txt
			echo "demulPath=$demulPath" >> $configPath/${recipeID}-${projID}_sampleID.txt
			echo "$ID" >> $configPath/${recipeID}-${projID}_sampleID.txt
			echo "### $configPath/${recipeID}-${projID}_sampleID.txt created ###"
		fi

		# Define delivery recipe and create delivery structure
		case $kit in
			Genome|Exome )
			dataType="DNA"
			echo "### dataType is $kit $dataType"
			for sampleID in `grep ID $configPath/${recipeID}-${projID}_sampleID.txt | cut -d= -f2`
			do
				mkdir -vp ${configPath}_RC1/$sampleID/{analysis/cna,bams,fastqs}
			done
			;;
			RNA )
			dataType="RNA"
			echo "### dataType is $kit $dataType"
			for sampleID in `grep ID $configPath/${recipeID}-${projID}_sampleID.txt | cut -d= -f2`
			do
				mkdir -vp ${configPath}_RC1/$sampleID/{analysis,bams,fastqs}
			done
			;;
			* )
			echo "### Kit type ($kit) not defined. Please check kitType."
		esac


		echo "### Running postrun delivery script on project ${recipeID}-${projID}"
		echo "sbatch --output=${configPath}/%x\_%j.out --job-name=${recipeID}-${projID}_postrunDelivery --export=RUNDIR=$runDir,RECIPEID=$recipeID,DEMULDIR=$demulPath,KITTYPE=$kitType,KIT=$kit,PROJID=$projID,CONFIGPATH=$configPath $slurmHome/postrun_Delivery.slurm"
		sbatch --output=${configPath}/%x\_%j.out --job-name=${recipeID}-${projID}_postrunDelivery --export=RUNDIR=$runDir,RECIPEID=$recipeID,DEMULDIR=$demulPath,KITTYPE=$kitType,KIT=$kit,PROJID=$projID,CONFIGPATH=$configPath $slurmHome/postrun_Delivery.slurm
		sleep 10
		echo "### Running postrun multiqc script on ${recipeID}-$projID"
		echo "sbatch --output=$runDir/multiqc/%x\_%j.out --job-name=${recipeID}-${projID}_multiqc --export=RUNDIR=$runDir,PROJID=$projID,RECIPEID=$recipeID $slurmHome/postrun_Multiqc.slurm"
		sbatch --output=$runDir/multiqc/%x\_%j.out --job-name=${recipeID}-${projID}_multiqc --export=RUNDIR=$runDir,PROJID=$projID,RECIPEID=$recipeID $slurmHome/postrun_Multiqc.slurm
	done
	sleep 10
done

echo "########################################################################"
echo "########################################################################"
echo "############################# End ######################################"
echo "########################################################################"
echo "########################################################################"
echo
echo
echo 
