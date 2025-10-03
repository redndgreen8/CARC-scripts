#!/bin/bash
#####################################################################
# Copyright (c) 2019 by University of Southern California,
# Institute of Translational Genomics; Authors: Yuxin Jin
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

##### Parameterized SLURM Script ####
#SBATCH --time=12:00:00
#SBATCH --mail-user=rdagnew@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --exclude=e20-41,e20-43,e20-45,e20-47
#SBATCH --partition=dtg 
#SBATCH --mem=0

module load usc python/3.9.2
module list

#Export necessary path
export PYTHONPATH=/project/davidwcr_264/Packages/pythonModule/python-3.9.2:$PYTHONPATH

#Define variables
configPath_DupRemoved=/home1/rdagnew/deliveryScripts/config/multiqc_config_KGP_DupRemoved.yaml
configPath=/home1/rdagnew/deliveryScripts/config/multiqc_config_KGP.yaml
multiqcPath=/project/davidwcr_264/Packages/pythonModule/python-3.9.2/bin/multiqc

projID=$PROJID
runDir=$RUNDIR
recipeID=$RECIPEID

echo "### projID=$projID ###"
echo "### runDir=$runDir ###"
echo "### recipeID=$recipeID ###"

#Command

echo "multiqcNormal=$multiqcPath -c $configPath -n ${recipeID}-$projID -o $runDir/multiqc $runDir/${recipeID}Out/$projID"
multiqcNormal="$multiqcPath -c $configPath -n ${recipeID}-${projID} -o $runDir/multiqc $runDir/${recipeID}Out/$projID"

echo "multiqcDupRemoved=$multiqcPath -c $configPath_DupRemoved -n ${recipeID}-${projID}_DupRemoved -o $runDir/multiqc $runDir/${recipeID}Out/$projID"
multiqcDupRemoved="$multiqcPath -c $configPath_DupRemoved -n ${recipeID}-${projID}_DupRemoved -o $runDir/multiqc $runDir/${recipeID}Out/$projID"


# Check output structure existance

if [[ ! -e ${runDir}/multiqc ]]; then
	echo "### multiqc folder doesn't exist"
	echo "### Start create multiqc folder"
	mkdir ${runDir}/multiqc
	echo "### Start run multiqc on DupRemoved $projID"
	$multiqcDupRemoved
	echo "### Start run multiqc on normal $projID"
	$multiqcNormal
else
	echo "### Found existing multiqc folder"
	echo "### Start looking for existing multiqc report"
	if [[ -f ${runDir}/multiqc/${recipeID}-${projID}.html ]] && [[ -f ${runDir}/multiqc/${recipeID}-${projID}_DupRemoved.html ]]; then
		echo "### multiqc might have already processed on $projID"
		echo " ### Please double check runDir"
		echo "###### Will exit!!! ######"
		exit
	elif [[ -f ${runDir}/multiqc/${recipeID}-${projID}.html ]] && [[ ! -f ${runDir}/multiqc/${recipeID}-${projID}_DupRemoved.html ]]; then
		echo "### Found multiqc report on normal $projID"
		echo "### Start run multiqc on DupRemoved $projID"
		$multiqcDupRemoved
	elif [[ ! -f ${runDir}/multiqc/${recipeID}-${projID}.html ]] && [[ -f ${runDir}/multiqc/${recipeID}-${projID}_DupRemoved.html ]]; then
		echo "### Found multiqc report on DupRemoved $projID"
		echo "### Start run multiqc on normal $projID"
		$multiqcNormal
	else
		echo "### No mutiqc report found"
		echo "### Start run multiqc on DupRemoved $projID"
        	$multiqcDupRemoved
		echo "### Start run multiqc on normal $projID"
		$multiqcNormal
	fi
fi
