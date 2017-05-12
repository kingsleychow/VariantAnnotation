#!/bin/bash
#PBS -r n

##Job settings
#PBS -N DCM_Reseq_BWA_GATK_run
#PBS -m be
#PBS -q bs_primary
#PBS -M b.jaydutt@nhcs.com.sg

##Job configuration
#PBS -V

##Job Resources
#PBS -l nodes=1:ppn=8
#PBS -l mem=30gb
####PBS -l walltime=120:00:00

annotations=/data1/seq_data/annotations/hg19
Ref=$annotations/genome.fa
dbSNP=$annotations/dbsnp_138.hg19.vcf
GATKs=/data1/apps/GATK/GenomeAnalysisTK-3.3-0
bwa=/data1/apps/bwa/bwa-0.7.10
bwa_index=$annotations/bwa_0_7_9a_index/genome.fa
picard=/data1/apps/picard/picard-tools-1.117
fastqc=/data1/apps/FastQC
samtools=/data1/apps/samtools/samtools-0.1.19

Java7=/usr/bin/java
memory=Xmx28g

bedtools=/data1/apps/bedtools/bedtools-2.17.0/bin

omni=$annotations/1000G_omni2.5.hg19.vcf
TenK=$annotations/1000G_phase1.snps.high_confidence.hg19.vcf
hapmap=$annotations/hapmap_3.3.hg19.vcf
mills=$annotations/Mills_and_1000G_gold_standard.indels.hg19.vcf
TenKIndel=$annotations/1000G_phase1.indels.hg19.vcf
snpEff=/data1/apps/snpEff/snpEff_4.1_L

tabix=/data1/apps/tabix-0.2.6

source /home/jaydutt/.bashrc

RunName=$1
Platform=$2
Lib_Kit=$3

data=/data1/seq_data/NHCS/$Platform/results

Coverage_Reports=$data/$RunName/Coverage_Report_v2

sftp=/data5/sftp/nhrisguest/incoming/upload/to_BRU

	mkdir $sftp/"$RunName"
	mkdir $sftp/$RunName/"Coverage_Report"


Coverage_script=$data/Scripts/$RunName/Coverage_Report


###########################################################################################

bash $Coverage_script/CoverageSummaryScript.sh $RunName $Platform $Lib_Kit

bash $Coverage_script/CallableCoverageMatrix.sh $RunName $Platform $Lib_Kit

###########################################################################################

## Coverage & Callable reports
	cp $Coverage_Reports/ProteinCodingTarget169gene/SummaryOutput_ProteinCodingTarget169gene_"$RunName".txt  $sftp/$RunName/Coverage_Report/SummaryOutput_ProteinCodingTarget_"$RunName".txt
	
	cp $Coverage_Reports/ProteinCodingTarget169gene/CallableByExon_ProteinCodingTarget169gene_"$RunName"_ICCNexteraV4.txt  $sftp/$RunName/Coverage_Report/CallableByExon_ProteinCodingTarget169gene_"$RunName"_ICCNexteraV4_169.txt

	cp $Coverage_Reports/ProteinCodingTarget169gene_30x/CallableByExon_ProteinCodingTarget169gene_30x_"$RunName"_ICCNexteraV4.txt  $sftp/$RunName/Coverage_Report/CallableByExon_ProteinCodingTarget169gene_30x_"$RunName"_ICCNexteraV4_169.txt

	cp $Coverage_Reports/CanonTranCodingTarget/CallableByExon_CanonTranCodingTarget_"$RunName"_ICCNexteraV4.txt  $sftp/$RunName/Coverage_Report/CallableByExon_CanonTranCodingTarget_"$RunName"_ICCNexteraV4_169.txt

	cp $Coverage_Reports/CanonTranCodingTarget_30x/CallableByExon_CanonTranCodingTarget_30x_"$RunName"_ICCNexteraV4.txt  $sftp/$RunName/Coverage_Report/CallableByExon_CanonTranCodingTarget_30x_"$RunName"_ICCNexteraV4_169.txt


