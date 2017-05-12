#!/bin/bash

######## Change RunName; FastqDir Location; assay & target values ##########

## RunName=150108_M02463_0029_000000000-AA4VF

RunName=$1
Platform=$2
Lib_Kit=$3

##From BaseSpace manually download
Data_dir=/data1/ngs/$Platform/$RunName
#Data_dir=/data1/seq_data/NHCS/NextSeq/$RunName
FastqDir=$Data_dir/Data/Intensities/BaseCalls

###From MiSeq machine
#FastqDir=/sequencer/MiSeq/Illumina/MiSeqOutput/$RunName/Data/Intensities/BaseCalls/

### Assay options: [Fluidigm | AgilentSureSelect |  AgilentHaloplex | Illumina]
assay=Illumina
### Target options: [Piezo |  PRDM16 | WARS2 | ICCv5 | ICCv6 | ICC_177_Dx | ICC_130_Dx | ICCNexteraV4]
target=$Lib_Kit
###########################################################################

SourceDir=/data1/seq_data/NHCS/$Platform/results
RunDir=/data1/seq_data/NHCS/$Platform/results/Reads
Scripts=/data1/seq_data/NHCS/$Platform/results/Scripts

#<<commentALL
mkdir $SourceDir/"$RunName"
mkdir $RunDir/"$RunName"
mkdir $Scripts/"$RunName"
mkdir $Scripts/$RunName/"BWA_GATK_Samtools"
mkdir $Scripts/$RunName/"Coverage_Report"


echo -e " samples"

### Copy the FastQ files into RunName ##ignore Undetermined, Negative Control, Positive Control Samples
SAMPLES=`ls $FastqDir | grep fastq.gz | awk 'BEGIN {FS="_S";} {print $1}'| sort -u| grep -v Undetermined | grep -v ^N`
echo -e " Samples copying from : $FastqDir "
for sample in $SAMPLES
do
	echo -e " Copied to : $RunDir/$RunName/"$sample""
	mkdir $RunDir/"$RunName"/"$sample"
	cp $FastqDir/"$sample"* $RunDir/"$RunName"/$sample/
	ls -ltr $RunDir/"$RunName"/$sample/
done

## Copy SampleSheet.csv
cp $Data_dir/SampleSheet.csv $SourceDir/$RunName/

if [ $assay == "Illumina" ] && [ $target == "ICCNexteraV4" ]
then
	cp -r /data1/seq_data/Store/Scripts/ICCNexteraV4/BWA_GATK_Samtools/*.sh $Scripts/$RunName/BWA_GATK_Samtools/
	cp -r /data1/seq_data/Store/Scripts/ICCNexteraV4/Coverage_Report/*.sh $Scripts/$RunName/Coverage_Report/
    sed -i '0,/RunName=/s/RunName=/RunName='$RunName'/' $Scripts/$RunName/BWA_GATK_Samtools/run*.sh
    sed -i '0,/Lib_Kit=/s/Lib_Kit=/Lib_Kit='$Lib_Kit'/' $Scripts/$RunName/BWA_GATK_Samtools/run*.sh
    sed -i '0,/Platform=/s/Platform=/Platform='$Platform'/' $Scripts/$RunName/BWA_GATK_Samtools/run*.sh
	echo -e "\n Assay=$assay \n Target=$target "
	echo -e "\nScripts copied into \n$Scripts/$RunName/BWA_GATK_Samtools/ \n$Scripts/$RunName/Coverage_Report/ \n "


elif [ $assay == "AgilentHaloplex" ] && [ $target == "MODY_panel" ]
then
	cp -r /data1/seq_data/NHCS/MiSeq/BWA_GATKv3.3_MODY_panel_pipeline_10Nov2014/BWA_GATK_Samtools/*.sh $Scripts/$RunName/BWA_GATK_Samtools/
 	cp -r /data1/seq_data/NHCS/MiSeq/BWA_GATKv3.3_MODY_panel_pipeline_10Nov2014/Coverage_Report/*.sh $Scripts/$RunName/Coverage_Report/
	sed -i '0,/RunName=/s/RunName=/RunName='$RunName'/' $Scripts/$RunName/BWA_GATK_Samtools/run*.sh
	sed -i '0,/RunName=/s/RunName=/RunName='$RunName'/' $Scripts/$RunName/BWA_GATK_Samtools/run*.sh
	sed -i '0,/Platform=/s/Platform=/Platform='$Platform'/' $Scripts/$RunName/BWA_GATK_Samtools/run*.sh
	echo -e "\n Assay=$assay \n Target=$target "
	echo -e "\nScripts copied into \n$Scripts/$RunName/BWA_GATK_Samtools/ \n$Scripts/$RunName/Coverage_Report/ \n "
	
else
	echo -e "\n!!!WARNING!!!\nAssay & Target names doesn't match"

fi



echo -e "\n  ##### Goto Scripts Dir & Run run***.sh script ######"
echo -e "\n cd "$Scripts"/"$RunName"/BWA_GATK_Samtools/"
echo -e "\n  ### Check the script files are ready to submit #####"
echo -e "\n  ##### Example bash run****.sh  ##### \n"



