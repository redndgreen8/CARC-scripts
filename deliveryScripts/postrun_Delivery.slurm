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
###SBATCH --output=%x\_%j.out
###SBATCH --job-name=postrunDelivery
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=dtg 
#SBATCH --mem=1GB

time=`date +%d-%m-%Y-%H-%M`

echo "########################################################################"
echo "Starting $0 at $time"

# Define variable
myhostname=`hostname`
runDir=$RUNDIR
recipeID=$RECIPEID
demulDir=$DEMULDIR
kitType=$KITTYPE
kit=$KIT
projID=$PROJID
configPath=$CONFIGPATH
releasePath=${CONFIGPATH}_RC1

echo "### myhostname=$myhostname ###"
echo "### runDir=$runDir ###"
echo "### recipeID=$recipeID ###"
echo "### demulDir=$demulDir ###"
echo "### kitType=$kitType ###"
echo "### kit=$kit ###"
echo "### projID=$projID ###"
echo "### configPath=$configPath ###"
echo "### releasePath=$releasePath ###"

echo "### --Running on $myhostname-- ###"

# Specify runfolder path
if [ -e $configPath/${recipeID}-${projID}_sampleID.txt ]; then
	echo "Found $configPath/${recipeID}-${projID}_sampleID.txt"
	echo "Start creating delivery structure."
else
	echo "### Missing $configPath/${recipeID}-${projID}_sampleID.txt"
	echo "### Please double check run directory path!!!"
	echo "### Exiting!!!"
	exit
fi


# Define variable

if [ -e $demulDir ]; then
	echo "Found $demulDir"
else
	echo "### Missing $demulDir"
	echo "### Please double check demulOut path"
	echo "### Exiting!!"
	exit
fi

# Define delivery recipe
case $kit in 
	Genome|Exome )
	dataType="DNA"
	cnas=exo.amp.tsv,exo.baf.png,exo.baf.seg.png,exo.baf.seg.vcf,exo.baf.tsv,exo.baf.txt,exo.cnaAbsBAF.png,exo.cnaBAF.png,exo.cna.seg,exo.cna.seg.vcf,exo.cna.tsv,exo.del.tsv,exo.hets.tsv,exo.png,exo.seg.png,exo_withhets.png,filt.amp.tsv,filt.baf.png,filt.baf.seg.png,filt.baf.txt,filt.cna.seg,filt.cna.seg.vcf,filt.cna.tsv,filt.del.tsv,filt.png,filt.seg.png,unfi.amp.tsv,unfi.baf.png,unfi.baf.seg.png,unfi.baf.txt,unfi.cna.seg,unfi.cna.seg.vcf,unfi.cna.tsv,unfi.del.tsv,unfi.png,unfi.seg.png
	analysis=.freebayes_All.snpEff.vcf,.HC_All.snpEff.vcf,.HC_All.vcf,.HC_All.vcf.idx,.manta.candidateSmallIndels.vcf.gz,.manta.candidateSmallIndels.vcf.gz.tbi,.manta.somaticSV.vcf.gz,.manta.somaticSV.vcf.gz.tbi,.REVseurat.snpEff.vcf,.REVseurat.vcf,.seurat.snpEff.vcf,.seurat.vcf,.strelka.passed.somatic.indels.snpEff.vcf,.strelka.passed.somatic.indels.vcf,.strelka.passed.somatic.snvs.snpEff.vcf,.strelka.passed.somatic.snvs.vcf,.trn.vcf,_Summary_allStats.txt,MuTect_All_Filtered.snpEff.vcf,MuTect_All_Filtered.vcf,MuTect_All_Filtered.vcf.idx
	bams=.bwa.final.bai,.bwa.final.bam,.bwa.final.bam.md5
	fastqs=.fastq.gz
	;;
	RNA )
	dataType="RNA"
	analysis=_genes.sf,_transcripts.sf,.rnaHC_All.snpEff.vcf,.rnaHC_All.vcf,.rnaHC_All.vcf.idx,.htSeqCounts,.junctions,.fusion_candidates.final.abridged,.fusion_predictions.abridged.tsv,.final-list_candidate-fusion-genes.caption.md.txt,.final-list_candidate-fusion-genes.GRCh37.txt,.final-list_candidate-fusion-genes.txt,.final-list_candidate-fusions.txt,.viruses_bacteria_phages.txt,.cufflinks.genes.fpkm_tracking,.cufflinks.isoforms.fpkm_tracking,.cufflinks.transcripts.gtf,.leafcutter.junc,.featureCounts.summary
	bams=.final.bai,.final.bam,.final.bam.md5
	fastqs=.fastq.gz
	;;
	* )
	echo "Kit type ($kit) not defined. Please check kitType."
esac

echo "Data type is $dataType"

# Start searching for files

if [ $dataType == "DNA" ]; then
	touch $configPath/cna.list
	touch $configPath/analysis.list
	touch $configPath/bams.list
	touch $configPath/fastqs.list
	
	for sample in `cat $configPath/${recipeID}-${projID}_sampleID.txt | grep ID | cut -d= -f2`
	do
		echo "################################################################"
		echo "### Start creating delivery structure for $dataType $sample ###"
		echo "################################################################"
		
		# Look for files deliver to cna folder
		for cna in `echo $cnas | tr "," "\n"`
		do
			echo "Start looking for $cna deliver to cna folder"
			find ${runDir}/${recipeID}Out/${projID} -type f -name "${sample}*${cna}" -exec echo {} \; >> $configPath/cnaCount.list
			cnaCount=`cat $configPath/cnaCount.list | wc -l`
			echo "Found $cnaCount $cna for $sample"
			cat $configPath/cnaCount.list >> $configPath/cna.list
			rm $configPath/cnaCount.list
		done
		echo "#################################################################"
		
		# Look for files deliver to analysis folder	
		for vcf in `echo $analysis | tr "," "\n"`
	do
			echo "Start looking for $vcf deliver to analysis folder"
			find ${runDir}/${recipeID}Out/${projID} -type f -name "${sample}*${vcf}" -exec echo {} \; >> $configPath/vcfCount.list
			vcfCount=`cat $configPath/vcfCount.list | wc -l`
			echo "Found $vcfCount $vcf for $sample"
			cat $configPath/vcfCount.list >> $configPath/analysis.list
			rm $configPath/vcfCount.list
		done
		echo "#################################################################"
		
		# Look for files deliver to bam folder
		for bam in `echo $bams | tr "," "\n"`
		do
			echo "Start looking for $bam deliver to bams folder"
			find ${runDir}/${recipeID}Out/${projID} -type f -name "${sample}*${bam}" -exec echo {} \; >> $configPath/bamCount.list
			bamCount=`cat $configPath/bamCount.list | wc -l`
			echo "Found $bamCount $bam for $sample"
			cat $configPath/bamCount.list >> $configPath/bams.list
			rm $configPath/bamCount.list
		done
		echo "#################################################################"
		
		# Look for files deliver to fastq folder
		for fastq in `echo $fastqs | tr "," "\n"`
		do
			echo "Start looking for $fastq deliver to fastqs folder"
			find $demulDir -type f -name "${sample}*${fastq}" -exec echo {} \; >> $configPath/fastqCount.list
			fastqCount=`cat $configPath/fastqCount.list | wc -l`
			echo "Found $fastqCount $fastq for $sample"
			cat $configPath/fastqCount.list >> $configPath/fastqs.list
			rm $configPath/fastqCount.list
		done
		echo "#################################################################"
	done
			
elif [ $dataType == "RNA" ]; then
	touch $configPath/analysis.list
	touch $configPath/bams.list
	touch $configPath/fastqs.list

	for sample in `cat $configPath/${recipeID}-${projID}_sampleID.txt | grep ID | cut -d= -f2`
	do
		echo "################################################################"
		echo "### Start creating delivery structure for $dataType $sample ###"
		echo "################################################################"


		# Look for files deliver to analysis folder
		for vcf in `echo $analysis | tr "," "\n"`
		do
			echo "Start looking for $vcf deliver to analysis folder"
			find ${runDir}/${recipeID}Out/${projID} -type f -name "${sample}*${vcf}" -exec echo {} \; >> $configPath/vcfCount.list
			vcfCount=`cat $configPath/vcfCount.list | wc -l`
			echo "Found $vcfCount $vcf for $sample"
			cat $configPath/vcfCount.list >> $configPath/analysis.list
			rm $configPath/vcfCount.list
		done
		echo "#################################################################"

		# Look for files deliver to bam folder
		for bam in `echo $bams | tr "," "\n"`
		do
			echo "Start looking for $bam deliver to bams folder"
			find ${runDir}/${recipeID}Out/${projID} -type f -name "${sample}*${bam}" -exec echo {} \; >> $configPath/bamCount.list
			bamCount=`cat $configPath/bamCount.list | wc -l`
			echo "Found $bamCount $bam for $sample"
			cat $configPath/bamCount.list >> $configPath/bams.list
			rm $configPath/bamCount.list
		done
		echo "#################################################################"

		# Look for files deliver to fastq folder
		for fastq in `echo $fastqs | tr "," "\n"`
		do
			echo "Start looking for $fastq deliver to fastqs folder"
			find $demulDir -type f -name "${sample}*${fastq}" -exec echo {} \; >> $configPath/fastqCount.list
			fastqCount=`cat $configPath/fastqCount.list | wc -l`
			echo "Found $fastqCount $fastq for $sample"
			cat $configPath/fastqCount.list >> $configPath/fastqs.list
			rm $configPath/fastqCount.list
		done
		echo "#################################################################"
	done
else
	echo "### Fail to define delivery recipe"
	echo "### Exiting!!!"
	exit
fi

# Create symlink for delivery
# Create bam link
echo "#################################################################"
echo "### Start creat links for delivery ###"
for bamFile in `cat $configPath/bams.list`
do
        targetDir=bams
        sampleID=`echo $bamFile | rev | cut -d/ -f1 | rev | cut -d_ -f1,2`
        bamTarget=`basename $bamFile`
        ln -sv $bamFile $releasePath/$sampleID/$targetDir/$bamTarget
done

sleep 10

# Create fastq link
for fastqFile in `cat $configPath/fastqs.list`
do
        targetDir=fastqs
        sampleID=`echo $fastqFile | rev | cut -d/ -f1 | rev | cut -d_ -f1,2`
        fastqTarget=`basename $fastqFile`
        ln -sv $fastqFile $releasePath/$sampleID/$targetDir/$fastqTarget
done

sleep 10

# Create analysis link
for analysisFile in `cat $configPath/analysis.list`
do
        targetDir=analysis
        sampleID=`echo $analysisFile | rev | cut -d/ -f1 | rev | cut -d_ -f1,2`
        analysisTarget=`basename $analysisFile`
        ln -sv $analysisFile $releasePath/$sampleID/$targetDir/$analysisTarget
done

sleep 10

if [ $dataType == "DNA" ]; then

# Create cna link
for cnaFile in `cat $configPath/cna.list`
do
	targetDir=analysis/cna
        sampleID=`echo $cnaFile | rev | cut -d/ -f1 | rev | cut -d_ -f1,2`
        cnaTarget=`basename $cnaFile`
        ln -sv $cnaFile $releasePath/$sampleID/$targetDir/$cnaTarget
done

sleep 10

elif [ $dataType == "RNA" ]; then

echo "No cna folder found for RNA"

fi



echo "Finishing $0 at $time"
echo "########################################################################"
echo "########################################################################"
echo "########################################################################"
echo "########################################################################"
echo "########################################################################"
