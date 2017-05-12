#!/bin/bash
#
#To modify script file by Sample 


###################################### Check these variables######################
##Run details
RunName=
Platform=
Lib_Kit=

ScriptFile=BWA_v0.7.12_GATK_v3.5_ICCNexteraV4_174Gene_v3.sh
##################################################################################

ReadFolder=Reads
ScriptFileLocation=/data1/seq_data/NHCS/$Platform/results/Scripts/$RunName/BWA_GATK_Samtools/

SampleDir=/data1/seq_data/NHCS/$Platform/results/$RunName
runDir=/data1/seq_data/NHCS/$Platform/results/$ReadFolder/$RunName
#Take out samples 'grep -v'
SAMPLES=`ls $runDir`
jobsToRun="$RunName"_"BWA_GATKJobs"
mkdir $ScriptFileLocation/"$jobsToRun"

#### Create a ScriptFile for each sample and SUBMIT the job into PBS nodes

#for sample in 10CD02790
for sample in $SAMPLES
do
	cp $ScriptFileLocation/$ScriptFile $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/Platform=/Platform='$Platform'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/ICC_/ICC_'$sample'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/for\ Lib\ in/for\ Lib\ in\ '$sample'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/RunName=/RunName='$RunName'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/Lib_Kit=/Lib_Kit='$Lib_Kit'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/stdout/'$sample'_stdout_'$Lib_Kit'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	sed -i 's/stderr/'$sample'_stderr_'$Lib_Kit'/' $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile
	cd $ScriptFileLocation/$jobsToRun/

	## Submit Jobs into PBS nodes
 	qsub -q normal $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile

	## To run the Job on master node
	#bash $ScriptFileLocation/$jobsToRun/"$sample"_$ScriptFile

	cd $ScriptFileLocation
	echo -e "\nSample :-  $sample submitted to Run"
done


