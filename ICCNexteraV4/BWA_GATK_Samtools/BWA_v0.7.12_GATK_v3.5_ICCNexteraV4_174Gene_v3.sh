#!/bin/bash
#PBS -r n

##Job settings
#PBS -N ICC_
#PBS -m be
#PBS -q bs_primary
#PBS -o stdout.$PBS_JOBID
#PBS -e stderr.$PBS_JOBID
#PBS -M b.jaydutt@nhcs.com.sg

##Job configuration
#PBS -V

##Job Resources
#PBS -l nodes=1:ppn=3
#PBS -l mem=15gb
#PBS -l walltime=120:00:00

################## These parameters will modify by the RUN script###############

RunName=
Platform=
Lib_Kit=

nt=4
dcovg=750
memory=Xmx16g

################## Check the target FIles are correct or NOT###################

##Target file from Ensembl70 & refSeq -- Overhangs exon and UTRs (+/-) 40bp or 300bp (for CNA analysis) into introns
onTarget_BED=/data1/seq_data/annotations/Target/ICCNexteraV4/ICCNexteraV3_CardiacDxPanel_174Genes.CDS.300bp_Padding.bed
onTarget_BED_169genes=/data1/seq_data/annotations/Target/ICCNexteraV4/ICCNexteraV2_CardiacDxPanel_169Genes.CDS.40bp_Padding.bed

##ProteinCodingregion
CDS_BED=/data1/seq_data/annotations/Target/ICCNexteraV4/ICCNexteraV3_CardiacDxPanel_174Genes.CDS.bed
CDS_BED_169genes=/data1/seq_data/annotations/Target/ICCNexteraV4/ICCNexteraV2_CardiacDxPanel_169Genes.CDS.bed
Canonical_Transcripts_169genes=/data1/seq_data/annotations/Target/ICCNexteraV4/ICCNexteraV3_169Genes_ProteinCoding_CanonicalTrans.mergeBed.bed

TopDir=/data1/seq_data/NHCS/$Platform/results
mkdir $TopDir/"$RunName"
mkdir $TopDir/$RunName/"$Lib_Kit"

SampleDir=/data1/seq_data/NHCS/$Platform/results/$RunName
runDir=/data1/seq_data/NHCS/$Platform/results/Reads/$RunName

sftp=/data5/sftp/nhrisguest/incoming/upload/to_BRU

	mkdir $sftp/"$RunName"
	mkdir $sftp/$RunName/"Coverage_Report"
	mkdir $sftp/$RunName/"HaplotypeCaller"
	mkdir $sftp/$RunName/"UnifiedGenotyper"

## create Output dir
outdir="gatk_snp_indel"
mkdir $SampleDir/$Lib_Kit/"$outdir"
MyOutDir=/data1/seq_data/NHCS/$Platform/results/$RunName/$Lib_Kit/$outdir

annotations=/data1/seq_data/annotations/hg19
Ref=$annotations/genome.fa
dbSNP=$annotations/dbsnp_138.hg19.vcf
GATKs=/data1/apps/GATK/GenomeAnalysisTK-3.5
bwa=/data1/apps/bwa/bwa-0.7.12
bwa_index=$annotations/bwa_0_7_12_index/genome.fa
picard=/data1/apps/picard/picard-tools-1.119
#fastqc=/data1/apps/FastQC/fastqc_v0.11.2
fastqc=/data1/apps/FastQC/fastqc_v0.11.4
samtools=/data1/apps/samtools/samtools-1.3
samtools=/data1/apps/samtools/bcftools-1.3
Java7=/usr/bin/java
bedtools=/data1/apps/bedtools/bedtools-2.17.0/bin
snpEff=/data1/apps/snpEff/snpEff_4.2

omni=$annotations/1000G_omni2.5.hg19.vcf
TenK=$annotations/1000G_phase1.snps.high_confidence.hg19.vcf
hapmap=$annotations/hapmap_3.3.hg19.vcf
mills=$annotations/Mills_and_1000G_gold_standard.indels.hg19.vcf
TenKIndel=$annotations/1000G_phase1.indels.hg19.vcf

tabix=/data1/apps/tabix-0.2.6

#source /home/jaydutt/.bashrc
. .bashrc
PATH=$PATH ; export PATH
################## These variables should modify based on runs###############

for Lib in
do
date

rm -r $runDir/$Lib/*prinseq*
rm -r $runDir/$Lib/*singleton*
rm -r $runDir/$Lib/*Trimed*
gunzip $runDir/$Lib/*.fastq.gz

SampleID=`ls $runDir/$Lib | grep R1_001.fastq | tail -1 | awk 'BEGIN {FS="_";} {print $1}'`
LaneID=`ls $runDir/$Lib | grep R1_001.fastq | tail -1 | awk 'BEGIN {FS="_";} {print $3}'`

cat $runDir/$Lib/*R1_001.fastq > $runDir/$Lib/"$SampleID"_R1.fastq
cat $runDir/$Lib/*R2_001.fastq > $runDir/$Lib/"$SampleID"_R2.fastq

Read1a=`ls $runDir/$Lib | grep R1.fastq | tail -1`
Read2a=`ls $runDir/$Lib | grep R2.fastq | tail -1`


perl /data1/apps/prinseq-lite-0.20.4/prinseq-lite.pl \
                -fastq $runDir/$Lib/$Read1a \
		-fastq2 $runDir/$Lib/$Read2a \
                -trim_qual_right 20 \
                -trim_qual_left 20 \
		-graph_data $runDir/$Lib/"$SampleID".gd\
		-graph_stats ld,gc,qd,de \
                -out_good $runDir/$Lib/"$SampleID"_Trimed


Read1=$runDir/$Lib/"$SampleID"_Trimed_1.fastq
Read2=$runDir/$Lib/"$SampleID"_Trimed_2.fastq

### Define Clean_BAM, onTarget_BAM, ProteinCoding_onTarget_BAM
onTarget_BAM="$SampleID".markDup.Realigned.recalibrated.OnTarget.q15.bam
Clean_BAM="$SampleID".markDup.Realigned.recalibrated.bam
ProteinCoding_onTarget_BAM="$SampleID".markDup.Realigned.recalibrated.ProteinCoding.OnTarget.q15.bam
ProteinCoding_onTarget_BAM_169genes="$SampleID".markDup.Realigned.recalibrated.169gene.ProteinCoding.OnTarget.q15.bam
Canonical_Transcripts_onTarget_BAM_169genes="$SampleID".markDup.Realigned.recalibrated.169gene.CanonTranCoding.OnTarget.q15.bam
Enrich_Report="Enrichment_Report"


### Create directory by SampleName/Target/ and SampleName/Exon and SampleName/Canonical_Transcripts
Target=Target
Exon=ProteinCodingTarget
Exon169gene=ProteinCodingTarget169gene
Canonical_Transcripts=CanonTranCodingTarget

mkdir $MyOutDir/"$SampleID"
mkdir $MyOutDir/"$SampleID"/"$Target"
mkdir $MyOutDir/"$SampleID"/"$Exon"
mkdir $MyOutDir/"$SampleID"/"$Exon169gene"
mkdir $MyOutDir/"$SampleID"/"$Canonical_Transcripts"

### Create directory  for EnrichmentReport and make symbolic link of files

mkdir $SampleDir/"$Enrich_Report"
mkdir $SampleDir/"$Enrich_Report"/"$SampleID"
mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$Target"
mkdir $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"

mkdir $sftp/$RunName/HaplotypeCaller/"$SampleID"
mkdir $sftp/$RunName/UnifiedGenotyper/"$SampleID"
cp $SampleDir/SampleSheet.csv $sftp/$RunName/SampleSheet.csv

### FastQC 
date

#<<FastQC
echo -e "FastQC started"
## $fastqc/fastqc $runDir/$Lib/$Read1
## $fastqc/fastqc $runDir/$Lib/$Read2

$fastqc/fastqc $Read1
$fastqc/fastqc $Read2

rm $runDir/$Lib/*prinseq*.fastq

bgzip $runDir/$Lib/"$SampleID"_R1.fastq
bgzip $runDir/$Lib/"$SampleID"_R2.fastq

unzip $runDir/$Lib/"$SampleID"_Trimed_1_fastqc.zip -d $runDir/$Lib
unzip $runDir/$Lib/"$SampleID"_Trimed_2_fastqc.zip -d $runDir/$Lib
mv $runDir/$Lib/*_1_fastqc.html $runDir/$Lib/"$SampleID"_Trimed_1_fastqc/
mv $runDir/$Lib/*_2_fastqc.html $runDir/$Lib/"$SampleID"_Trimed_2_fastqc/

#FastQC

$bwa/bwa 2> /data1/tmp/bwa.txt
bwaver=`cat /data1/tmp/bwa.txt | grep Version | awk '{print $2}'`

### BWA alignment
### -M	 Mark shorter split hits as secondary (for Picard compatibility). -t number of threads
echo -e "BWA mapping started"
## $bwa/bwa mem -M -t $nt $bwa_index $runDir/$Lib/$Read1 $runDir/$Lib/$Read2 > $MyOutDir/"$SampleID"/$SampleID.sam
$bwa/bwa mem -M -t $nt $bwa_index $Read1 $Read2 > $MyOutDir/"$SampleID"/$SampleID.sam

## SAM to BAM file
samtools view -bS $MyOutDir/"$SampleID"/"$SampleID".sam > $MyOutDir/"$SampleID"/"$SampleID".bam1
samtools view -h -F 4 -b $MyOutDir/"$SampleID"/"$SampleID".bam1 > $MyOutDir/"$SampleID"/"$SampleID".bam

samtools flagstat $MyOutDir/"$SampleID"/"$SampleID".bam1 > $MyOutDir/"$SampleID"/"$SampleID".TotalReads.flagstat
samtools flagstat $MyOutDir/"$SampleID"/"$SampleID".bam > $MyOutDir/"$SampleID"/"$SampleID".mapped.flagstat

##Remove SAM file
rm -rf $MyOutDir/"$SampleID"/"$SampleID".sam
rm -rf $MyOutDir/"$SampleID"/"$SampleID".bam1

date
echo -e "BWA run completed"

bgzip $runDir/$Lib/"$SampleID"_Trimed_1.fastq
bgzip $runDir/$Lib/"$SampleID"_Trimed_2.fastq

### Add ReadGroup tag
java -$memory -jar -Djava.io.tmpdir=/data1/tmp/ $picard/AddOrReplaceReadGroups.jar \
	INPUT=$MyOutDir/"$SampleID"/"$SampleID".bam \
	OUTPUT=$MyOutDir/"$SampleID"/"$SampleID".bam.bam \
	SORT_ORDER=coordinate \
	RGID=$RunName \
	RGLB=$SampleID"_150PEBC" \
	RGPL=ILLUMINA \
	RGPU=$LaneID \
	RGSM=$SampleID \
	RGCN=NHRIS \
	RGDS=BWA-MEM."$bwaver" \
	VALIDATION_STRINGENCY=SILENT

##  Allowable options for 'RGPL'are ILLUMINA,SLX,SOLEXA,SOLID,454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN

mv $MyOutDir/"$SampleID"/"$SampleID".bam.bam $MyOutDir/"$SampleID"/"$SampleID".bam
samtools index $MyOutDir/"$SampleID"/"$SampleID".bam


################ Marking duplicates.
java -Djava.io.tmpdir=/data1/tmp/ -$memory -jar $picard/MarkDuplicates.jar \
	INPUT=$MyOutDir/"$SampleID"/"$SampleID".bam \
	OUTPUT=$MyOutDir/"$SampleID"/"$SampleID".markDup.bam \
	METRICS_FILE=$MyOutDir/"$SampleID"/"$SampleID".markDup.bam.metrics.txt \
	VALIDATION_STRINGENCY=LENIENT

samtools flagstat $MyOutDir/"$SampleID"/"$SampleID".markDup.bam > $MyOutDir/"$SampleID"/"$SampleID".MarkDupReads.flagstat
		
samtools index $MyOutDir/"$SampleID"/"$SampleID".markDup.bam

###Creating Intervals 
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-nt $nt \
	-I $MyOutDir/"$SampleID"/"$SampleID".markDup.bam \
	-R $Ref \
	-L $onTarget_BED \
	--out $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals \
        -known $mills \
	-known $TenKIndel
	echo -e "\nDone Creating Intervals : $Lib"


###Indel Realigning  ####################### #####  
$Java7 -Djava.io.tmpdir=/data1/tmp/ -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-I $MyOutDir/"$SampleID"/"$SampleID".markDup.bam \
	-R $Ref \
	-T IndelRealigner \
	-targetIntervals $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals \
	--out $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam \
	-known $mills \
	-known $TenKIndel \
	-model USE_READS

echo -e "\nDone Indel Realigning : $Lib"
samtools index $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam


###CountCovariates
$Java7 -Djava.io.tmpdir=/data1/tmp/ -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-I $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam \
	--knownSites $dbSNP \
	--knownSites $mills \
	--knownSites $TenKIndel \
	-T BaseRecalibrator \
	-nct $nt \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov ContextCovariate \
	--out $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam.recalibTable.grp
echo -e "\nDone : CountCovariates : $Lib"

##TableRecalibration 
$Java7 -Djava.io.tmpdir=/data1/tmp/ -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-I $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam \
	-T PrintReads \
	-nct $nt \
	--out $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam \
	-BQSR $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam.recalibTable.grp \
	-S LENIENT

cp $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bai $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam.bai
####################################################################################

### Generate read on $Target MAP qual > 15
samtools view -uq 15 $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam | $bedtools/intersectBed -abam stdin -b $onTarget_BED -u > $MyOutDir/"$SampleID"/$Target/"$SampleID".markDup.Realigned.recalibrated.OnTarget.q15.bam
samtools index $MyOutDir/"$SampleID"/$Target/"$SampleID".markDup.Realigned.recalibrated.OnTarget.q15.bam

### Generate read on Protein coding target  MAP qual > 15
samtools view -uq 15 $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam | $bedtools/intersectBed -abam stdin -b $CDS_BED -u > $MyOutDir/"$SampleID"/$Exon/"$SampleID".markDup.Realigned.recalibrated.ProteinCoding.OnTarget.q15.bam
samtools index $MyOutDir/"$SampleID"/$Exon/"$SampleID".markDup.Realigned.recalibrated.ProteinCoding.OnTarget.q15.bam

### Generate read on Protein coding based on 169 gene CDS target  MAP qual > 15
samtools view -uq 15 $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam | $bedtools/intersectBed -abam stdin -b $CDS_BED_169genes -u > $MyOutDir/"$SampleID"/$Exon169gene/"$SampleID".markDup.Realigned.recalibrated.169gene.ProteinCoding.OnTarget.q15.bam
samtools index $MyOutDir/"$SampleID"/$Exon169gene/"$SampleID".markDup.Realigned.recalibrated.169gene.ProteinCoding.OnTarget.q15.bam

### Generate read on Protein coding based on 169 Canonical Transcripts CDS target  MAP qual > 15
samtools view -uq 15 $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam | $bedtools/intersectBed -abam stdin -b $Canonical_Transcripts_169genes -u > $MyOutDir/"$SampleID"/$Canonical_Transcripts/"$SampleID".markDup.Realigned.recalibrated.169gene.CanonTranCoding.OnTarget.q15.bam
samtools index $MyOutDir/"$SampleID"/$Canonical_Transcripts/"$SampleID".markDup.Realigned.recalibrated.169gene.CanonTranCoding.OnTarget.q15.bam

### Generating flagstat reports for ###Take out chimeric reads to get exact number of reads in flagstat -F 0x100
samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Target/$onTarget_BAM | samtools flagstat - > $MyOutDir/"$SampleID"/$Target/$onTarget_BAM.flagstat
samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM | samtools flagstat - > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.flagstat
samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes | samtools flagstat - > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.flagstat
samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes | samtools flagstat - > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.flagstat

samtools view -F 0x100 -uq 15 $MyOutDir/"$SampleID"/$Clean_BAM | samtools flagstat - > $MyOutDir/"$SampleID"/$Clean_BAM.q15.flagstat
samtools view -u -F 0x100 $MyOutDir/"$SampleID"/$Clean_BAM | samtools flagstat - > $MyOutDir/"$SampleID"/$Clean_BAM.flagstat

####################################################################################
####################################################################################
CDS_BED

### CoverageStas [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM | $bedtools/coverageBed -abam stdin -b $CDS_BED  > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats
### CoveragePerBase [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM | $bedtools/coverageBed -abam stdin -b $CDS_BED -d > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoveragePerBase
### CoverageHistogram [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM | $bedtools/coverageBed -abam stdin -b $CDS_BED -hist > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageHist

####################################################################################
CDS_BED_169genes

### CoverageStas [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes | $bedtools/coverageBed -abam stdin -b $CDS_BED_169genes  > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats
### CoveragePerBase [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes | $bedtools/coverageBed -abam stdin -b $CDS_BED_169genes -d > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoveragePerBase
### CoverageHistogram [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes | $bedtools/coverageBed -abam stdin -b $CDS_BED_169genes -hist > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageHist

####################################################################################
Canonical_Transcripts_169genes

samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes | $bedtools/coverageBed -abam stdin -b $Canonical_Transcripts_169genes  > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats
### CoveragePerBase [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes | $bedtools/coverageBed -abam stdin -b $Canonical_Transcripts_169genes -d > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoveragePerBase
### CoverageHistogram [ duplicateReads were removed to calculate stats]
samtools view -uF 0x400 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes | $bedtools/coverageBed -abam stdin -b $Canonical_Transcripts_169genes -hist > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageHist

####################################################################################
####################################################################################


CDS_BED
### sort the CoverageBed, CoverageStats output by chr,position
sort -k 1,1 -k 2,2n $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats.sorted
mv $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats.sorted $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats

echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoverageStats completed : $Lib"

### sort the CoverageBed,CoveragePerBase output by chr,position,coverage
sort -k 1,1 -k 2,2n -k 5,5n $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoveragePerBase  > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoveragePerBase.sorted
mv $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoveragePerBase.sorted $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoveragePerBase
echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoveragePerBase completed : $Lib"

####################################################################################

CDS_BED_169genes
### sort the CoverageBed, CoverageStats output by chr,position
sort -k 1,1 -k 2,2n $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats.sorted
mv $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats.sorted $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats

echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoverageStats completed : $Lib"

### sort the CoverageBed,CoveragePerBase output by chr,position,coverage
sort -k 1,1 -k 2,2n -k 5,5n $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoveragePerBase  > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoveragePerBase.sorted
mv $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoveragePerBase.sorted $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoveragePerBase
echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoveragePerBase completed : $Lib"
####################################################################################

Canonical_Transcripts_169genes
### sort the CoverageBed, CoverageStats output by chr,position
sort -k 1,1 -k 2,2n $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats.sorted
mv $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats.sorted $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats

echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoverageStats completed : $Lib"

### sort the CoverageBed,CoveragePerBase output by chr,position,coverage
sort -k 1,1 -k 2,2n -k 5,5n $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoveragePerBase  > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoveragePerBase.sorted
mv $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoveragePerBase.sorted $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoveragePerBase
echo -e "\nDone : Sorting of Protein coding OnTarget (Targets).CoveragePerBase completed : $Lib"

####################################################################################
####################################################################################

CDS_BED
#################################
###  Number of callable Bases by CDS at 10x
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM \
	-o $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 10 \
	-l INFO -format BED \
	-L $CDS_BED \
	-summary $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"

###  Number of callable Bases by CDS at 20x

mkdir $MyOutDir/"$SampleID"/"ProteinCodingTarget_20x"
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM \
	-o $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 20 \
	-l INFO -format BED \
	-L $CDS_BED \
	-summary $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"

###  Number of callable Bases by CDS at 30x

mkdir $MyOutDir/"$SampleID"/"ProteinCodingTarget_30x"
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM \
	-o $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 30 \
	-l INFO -format BED \
	-L $CDS_BED \
	-summary $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"


### Depth of coverage in CDS
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T DepthOfCoverage \
	-I $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM \
	-o $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.DepthOfCvg.txt  \
	-L $CDS_BED \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--includeDeletions -ct 5 -ct 10 -ct 20 -ct 30 -ct 50 -ct 100\
	--outputFormat table --omitDepthOutputAtEachBase

echo -e " \n Done : DepthOfCoverage for  : $Lib"

### Getting mean coverage per Exon
	cat $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.DepthOfCvg.txt.sample_interval_summary | sed 's/chr16:997401/chr16:997401-997401/' | awk 'BEGIN {FS=":"} {OFS="\t"} {print $1,$2}'| awk 'BEGIN {FS="-"} {OFS="\t"} {print $1,$2}' |  awk 'BEGIN {FS=" "} {OFS="\t"} {print $1,$2-1,$3,$4,$5,$9}' | sed '1d'| $bedtools/sortBed -i stdin > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.temp.depth1.bed
	cut -f4 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats | awk 'BEGIN {FS="|";} {print $1}' > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.temp.depth2.bed 
	paste <(cut -f1-3 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.temp.depth1.bed) <(cut -f4 $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.temp.depth2.bed ) <(cut -f4- $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.temp.depth1.bed) | cut -f1-4,6 > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed 

	rm -rf $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.temp*
	ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $MyOutDir/"$SampleID"/"$Exon"_20x/
	ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $MyOutDir/"$SampleID"/"$Exon"_30x/

####################################################################################

CDS_BED_169genes
#################################
###  Number of callable Bases by CDS_169 at 10x
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 10 \
	-l INFO -format BED \
	-L $CDS_BED_169genes \
	-summary $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"

###  Number of callable Bases by CDS_169 at 20x

mkdir $MyOutDir/"$SampleID"/"ProteinCodingTarget169gene_20x"
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 20 \
	-l INFO -format BED \
	-L $CDS_BED_169genes \
	-summary $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"

###  Number of callable Bases by CDS_169 at 30x

mkdir $MyOutDir/"$SampleID"/"ProteinCodingTarget169gene_30x"
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 30 \
	-l INFO -format BED \
	-L $CDS_BED_169genes \
	-summary $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"


### Depth of coverage in CDS_169
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T DepthOfCoverage \
	-I $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.DepthOfCvg.txt  \
	-L $CDS_BED_169genes \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--includeDeletions -ct 5 -ct 10 -ct 20 -ct 30 -ct 50 -ct 100\
	--outputFormat table --omitDepthOutputAtEachBase

echo -e " \n Done : DepthOfCoverage for  : $Lib"

### Getting mean coverage per Exon
	cat $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary | sed 's/chr16:997401/chr16:997401-997401/' | awk 'BEGIN {FS=":"} {OFS="\t"} {print $1,$2}'| awk 'BEGIN {FS="-"} {OFS="\t"} {print $1,$2}' |  awk 'BEGIN {FS=" "} {OFS="\t"} {print $1,$2-1,$3,$4,$5,$9}' | sed '1d'| $bedtools/sortBed -i stdin > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.temp.depth1.bed
	cut -f4 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats | awk 'BEGIN {FS="|";} {print $1}' > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.temp.depth2.bed 
	paste <(cut -f1-3 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.temp.depth1.bed) <(cut -f4 $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.temp.depth2.bed ) <(cut -f4- $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.temp.depth1.bed) | cut -f1-4,6 > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed 

	rm -rf $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.temp*
	ln -s $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $MyOutDir/"$SampleID"/"$Exon169gene"_20x/
	ln -s $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $MyOutDir/"$SampleID"/"$Exon169gene"_30x/

####################################################################################

Canonical_Transcripts
#################################
###  Number of callable Bases by Canonical_Transcripts at 10x
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 10 \
	-l INFO -format BED \
	-L $Canonical_Transcripts_169genes \
	-summary $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"

###  Number of callable Bases by CDS_169 at 20x

mkdir $MyOutDir/"$SampleID"/"CanonTranCodingTarget_20x"
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 20 \
	-l INFO -format BED \
	-L $Canonical_Transcripts_169genes \
	-summary $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"

###  Number of callable Bases by CDS_169 at 30x

mkdir $MyOutDir/"$SampleID"/"CanonTranCodingTarget_30x"
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T CallableLoci \
	-I $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--minDepth 30 \
	-l INFO -format BED \
	-L $Canonical_Transcripts_169genes \
	-summary $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.summary
echo -e " \n Done : CallableLoci  for : $Lib"


### Depth of coverage in CDS
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T DepthOfCoverage \
	-I $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes \
	-o $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.DepthOfCvg.txt  \
	-L $Canonical_Transcripts_169genes \
	--minMappingQuality 20 \
	--minBaseQuality 20 \
	--includeDeletions -ct 5 -ct 10 -ct 20 -ct 30 -ct 50 -ct 100\
	--outputFormat table --omitDepthOutputAtEachBase

echo -e " \n Done : DepthOfCoverage for  : $Lib"

### Getting mean coverage per Exon
	cat $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary | sed 's/chr16:997401/chr16:997401-997401/' | awk 'BEGIN {FS=":"} {OFS="\t"} {print $1,$2}'| awk 'BEGIN {FS="-"} {OFS="\t"} {print $1,$2}' |  awk 'BEGIN {FS=" "} {OFS="\t"} {print $1,$2-1,$3,$4,$5,$9}' | sed '1d'| $bedtools/sortBed -i stdin > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.temp.depth1.bed
	cut -f4 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats | awk 'BEGIN {FS="|";} {print $1}' > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.temp.depth2.bed 
	paste <(cut -f1-3 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.temp.depth1.bed) <(cut -f4 $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.temp.depth2.bed ) <(cut -f4- $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.temp.depth1.bed) | cut -f1-4,6 > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed 

	rm -rf $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.temp*
	ln -s $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/
	ln -s $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.DepthOfCvg.txt.sample_interval_summary.MeanCvg.bed $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/

####################################################################################
####################################################################################

### CollectAlignmentSummaryMetrics by CDS
java -$memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
	INPUT=$MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM \
	OUTPUT=$MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.AlignSumMet.txt \
	REFERENCE_SEQUENCE=$Ref \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT
echo -e " \n Done : CollectAlignmentSummaryMetrics for  : $Lib"

### CollectAlignmentSummaryMetrics by Target
java -$memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
	INPUT=$MyOutDir/"$SampleID"/$Target/$onTarget_BAM \
	OUTPUT=$MyOutDir/"$SampleID"/$Target/$onTarget_BAM.AlignSumMet.txt \
	REFERENCE_SEQUENCE=$Ref \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT
echo -e " \n Done : CollectAlignmentSummaryMetrics for  : $Lib"

### CollectAlignmentSummaryMetrics by CDS_169
java -$memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
	INPUT=$MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes \
	OUTPUT=$MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.AlignSumMet.txt \
	REFERENCE_SEQUENCE=$Ref \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT
echo -e " \n Done : CollectAlignmentSummaryMetrics for  : $Lib"

### CollectAlignmentSummaryMetrics by CDS Canonical Transcripts
java -$memory -jar $picard/CollectAlignmentSummaryMetrics.jar \
	INPUT=$MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes \
	OUTPUT=$MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.AlignSumMet.txt \
	REFERENCE_SEQUENCE=$Ref \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT
echo -e " \n Done : CollectAlignmentSummaryMetrics for  : $Lib"
####################################################################################
####################################################################################

CDS_BED
### Generate Bases callable by CDS_BED
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable | intersectBed -a $CDS_BED -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.byTarget

### Generate Bases callable by CDS_BED at 20x
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable | intersectBed -a $CDS_BED -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget > $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/"$Exon"_20x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget

### Generate Bases callable by CDS_BED at 30x
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable | intersectBed -a $CDS_BED -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget > $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/"$Exon"_30x/$ProteinCoding_onTarget_BAM.bases.callable.byTarget

####################################################################################

CDS_BED_169genes
### Generate Bases callable by CDS_BED_169genes
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable | intersectBed -a $CDS_BED_169genes -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget

### Generate Bases callable by CDS_BED at 20x
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable | intersectBed -a $CDS_BED_169genes -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget > $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/"$Exon169gene"_20x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget

### Generate Bases callable by CDS_BED at 30x
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable | intersectBed -a $CDS_BED_169genes -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget > $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/"$Exon169gene"_30x/$ProteinCoding_onTarget_BAM_169genes.bases.callable.byTarget

####################################################################################

Canonical_Transcripts
### Generate Bases callable by Canonical_Transcripts
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable | intersectBed -a $Canonical_Transcripts_169genes -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget > $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget

### Generate Bases callable by CDS_BED at 20x
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable | intersectBed -a $Canonical_Transcripts_169genes -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget > $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget

### Generate Bases callable by CDS_BED at 30x
#http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html
awk 'BEGIN{FS=" "; OFS="\t";} {print $1,$2,$3,$4}' $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable | intersectBed -a $Canonical_Transcripts_169genes -b stdin -wao | awk '{OFS="\t"; print $1,$2,$3,$4,$8"="$9}' | /data1/apps/filo-master/bin/groupBy -i stdin -grp 1,2,3,4 -c 5 -ops collapse > $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget

#### Perl script will generate this format 
##Chr     Start   End     Gene    CALLABLE        LOW_COVERAGE    NO_COVERAGE
##chr1    237205822       237205869       RYR2|1|ENSG00000198626  47      0       0
perl /data1/seq_data/Store/Scripts/CallableLoci_rearrage_Coverage.pl $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget > $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget.temp
mv $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget.temp $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/$Canonical_Transcripts_onTarget_BAM_169genes.bases.callable.byTarget

####################################################################################
#####################################################################################

ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats $MyOutDir/"$SampleID"/"$Exon"_20x/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats $MyOutDir/"$SampleID"/"$Exon"_30x/
ln -s $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats $MyOutDir/"$SampleID"/"$Exon169gene"_20x/
ln -s $MyOutDir/"$SampleID"/$Exon169gene/$ProteinCoding_onTarget_BAM_169genes.CoverageStats $MyOutDir/"$SampleID"/"$Exon169gene"_30x/
ln -s $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_20x/
ln -s $MyOutDir/"$SampleID"/$Canonical_Transcripts/$Canonical_Transcripts_onTarget_BAM_169genes.CoverageStats $MyOutDir/"$SampleID"/"$Canonical_Transcripts"_30x/

####################################################################################

### Generate  soft links in /data/results/$Platform/$RunName/$Enrich_Report folder

ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.byTarget $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Clean_BAM.q15.flagstat $SampleDir/"$Enrich_Report"/$SampleID/
ln -s $MyOutDir/"$SampleID"/$Clean_BAM.flagstat $SampleDir/"$Enrich_Report"/$SampleID/

ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.flagstat $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageStats $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoveragePerBase $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.CoverageHist $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.bases.callable.summary $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.AlignSumMet.txt $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/
ln -s $MyOutDir/"$SampleID"/$Exon/$ProteinCoding_onTarget_BAM.DepthOfCvg.txt.sample_interval_summary $SampleDir/"$Enrich_Report"/$SampleID/"$Exon"/

#######################################################################################################################################
# ************************************************GATK - Indels and SNPs calling************************************************************#
########################################################################################################################################
mkdir $MyOutDir/$SampleID/$Target/"UnifiedGenotyper"

###Calling UnifiedGenotyper  ###
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar -l INFO \
	-R $Ref \
	-L $onTarget_BED \
	-I $MyOutDir/"$SampleID"/$Target/$onTarget_BAM \
	-T UnifiedGenotyper \
	-baq CALCULATE_AS_NECESSARY \
	-minIndelCnt 4 \
	-glm BOTH \
	--dbsnp $dbSNP \
	-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.indel.vcf \
	-A Coverage \
	-A AlleleBalance \
	-G Standard \
	-gt_mode DISCOVERY \
	--min_base_quality_score 20 \
	--metrics_file $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.indel.metrics \
	--num_threads $nt \
	-stand_call_conf 30.0 \
	-stand_emit_conf 10.0 \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	--computeSLOD

echo -e "\nDone : SNP & Indels calling : $Lib"

#### print Indels
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-T SelectVariants \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-R $Ref \
	-nt $nt \
	--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.indel.vcf \
	--selectTypeToInclude INDEL \
	-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.indel.vcf


echo -e "\nDone : Indels calling : $Lib"


### Filter Indels
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantFiltration \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.indel.vcf \
	-baq CALCULATE_AS_NECESSARY \
	--genotypeFilterExpression "GQ < 20.0" \
	--filterExpression "QD < 2.0" \
	--filterExpression "ReadPosRankSum < -20.0" \
	--filterExpression "FS > 200.0" \
	--genotypeFilterName lowGQ \
	--filterName QDFilter \
	--filterName ReadPosFilter \
	--filterName FSFilter

#### Print SNVs
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-T SelectVariants \
 	-R $Ref \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-nt $nt \
 	--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.indel.vcf \
	--selectTypeToInclude SNP \
 	-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.vcf


echo -e "\nDone : Indel filtering : $Lib"
### Filter SNVs 
$Java7 -$memory -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantFiltration \
	--downsample_to_coverage $dcovg \
	--downsampling_type BY_SAMPLE \
	-o $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.vcf \
	--mask $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.indel.vcf \
	--maskName InDel \
	--genotypeFilterExpression "GQ < 20.0" \
	--filterExpression "QD < 2.0" \
	--filterExpression "MQ < 40.0" \
	--filterExpression "FS > 60.0" \
	--filterExpression "MQRankSum < -12.5" \
	--filterExpression "ReadPosRankSum < -8.0" \
	--genotypeFilterName lowGQ \
	--filterName QDFilter \
	--filterName MQFilter \
	--filterName FSFilter \
	--filterName MQRankSumFilter \
	--filterName ReadPosFilter

echo -e "\nDone : SNP filtering : $Lib"

### Filter SNV's variant evaluation 
$Java7 -$memory -jar $snpEff/SnpSift.jar tstv $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf > $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf.all.snp.TsTv.txt

ln -s $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf.all.snp.TsTv.txt $SampleDir/"$Enrich_Report"/$SampleID/

echo -e "\nDone :  Filter SNV's variant evaluation : $Lib"

echo -e "\nDone : GATK & SAMTools analysis were completed for : $Lib"
###Generating final filtered SNP VCF file
ln -s $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.final.UnifiedGenotyper.snp.vcf
###Generating final filtered Indel VCF file
ln -s $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.final.UnifiedGenotyper.indel.vcf


####  HaplotypeCaller  #############################
############## GATK new variant Call method -T HaplotypeCaller; This will replace UnifiedGenotyper in future############
mkdir $MyOutDir/"$SampleID"/$Target/"HaplotypeCaller"
## HaplotypeCaller
$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-L $onTarget_BED \
	-I $MyOutDir/"$SampleID"/$Target/$onTarget_BAM \
	-T HaplotypeCaller \
	-A Coverage \
	-A AlleleBalance \
	-G Standard \
	-gt_mode DISCOVERY \
	-nct $nt \
	--dbsnp $dbSNP \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.vcf \
	-stand_call_conf 30 \
	-stand_emit_conf 10

##
#### Variant Annotater for haplotypescore Allele Balance
$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantAnnotator \
	-L $onTarget_BED \
	-I $MyOutDir/"$SampleID"/$Target/$onTarget_BAM \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.varAnnotate.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.vcf \
	-A Coverage \
	-A AlleleBalance \
	-A HaplotypeScore \
	-A InbreedingCoeff \
	-A HomopolymerRun \
	-A HardyWeinberg \
	-A GCContent \
	--dbsnp $dbSNP \
	-nt $nt 


#### Print Indels
$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $Ref \
	-nt $nt \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.varAnnotate.vcf \
	--selectTypeToInclude INDEL \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.indel.vcf

echo -e "\nDone : HaplotypeCaller Indels calling : $Lib"


### Filter Indels
$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantFiltration \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.indel.vcf \
	-baq CALCULATE_AS_NECESSARY \
	--genotypeFilterExpression "GQ < 20.0" \
	--filterExpression "QD < 2.0" \
	--filterExpression "ReadPosRankSum < -20.0" \
	--filterExpression "FS > 200.0" \
	--genotypeFilterName lowGQ \
	--filterName QDFilter \
	--filterName ReadPosFilter \
	--filterName FSFilter

#### Print SNVs
$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-T SelectVariants \
 	-R $Ref \
	-nt $nt \
 	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.varAnnotate.vcf \
	--selectTypeToInclude SNP \
 	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.vcf

echo -e "\nDone : HaplotypeCaller Indel filtering : $Lib"

### Filter SNVs 
$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-T VariantFiltration \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf \
	--variant $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.vcf \
	--mask $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.indel.vcf \
	--maskName InDel \
	--genotypeFilterExpression "GQ < 20.0" \
	--filterExpression "QD < 2.0" \
	--filterExpression "MQ < 40.0" \
	--filterExpression "FS > 60.0" \
	--filterExpression "MQRankSum < -12.5" \
	--filterExpression "ReadPosRankSum < -8.0" \
	--genotypeFilterName lowGQ \
	--filterName QDFilter \
	--filterName MQFilter \
	--filterName FSFilter \
	--filterName MQRankSumFilter \
	--filterName ReadPosFilter

echo -e "\nDone :  HaplotypeCaller SNP filtering : $Lib"

### Filter SNV's variant evaluation 
$Java7 -$memory -jar $snpEff/SnpSift.jar tstv $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf > $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf.all.snp.TsTv.txt

ln -s $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf.all.snp.TsTv.txt $SampleDir/"$Enrich_Report"/$SampleID/

echo -e "\nDone :  Filter SNV's variant evaluation : $Lib"
echo -e "\nDone : GATK & SAMTools analysis were completed for : $Lib"

###Generating final filtered SNP VCF file
ln -s $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.final.HaplotypeCaller.snp.vcf
###Generating final filtered Indel VCF file
ln -s $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.final.HaplotypeCaller.indel.vcf

####  HaplotypeCaller  #############################
############## GATK new variant Call method -T HaplotypeCaller; This will replace UnifiedGenotyper in future############
## HaplotypeCaller
## HaplotypeCaller in gVCF mode

$Java7 -Xmx8g -jar $GATKs/GenomeAnalysisTK.jar \
	-R $Ref \
	-L $onTarget_BED \
	-I $MyOutDir/"$SampleID"/$Target/$onTarget_BAM \
	-T HaplotypeCaller \
	--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
	-A Coverage \
	-A AlleleBalance \
	-G Standard \
	-nct $nt \
	--dbsnp $dbSNP \
	-o $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf
	-stand_call_conf 30 \
	-stand_emit_conf 10


####  Data Transfer to SFTP   #############################

	mkdir $sftp/$RunName/HaplotypeCaller/"$SampleID"
	mkdir $sftp/$RunName/UnifiedGenotyper/"$SampleID"

	cp $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf  $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf

	cp $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf  $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf

	cp $MyOutDir/"$SampleID"/$Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf  $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf

##VCF files UnifiedGenotyper
	cp $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf  $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf

	cp $MyOutDir/"$SampleID"/$Target/UnifiedGenotyper/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf  $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf

	bgzip $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf
	bgzip $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf
	bgzip $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf
	bgzip $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf
	bgzip $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf

	$tabix/tabix -p vcf $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf.gz
	$tabix/tabix -p vcf $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf.gz
	$tabix/tabix -p vcf $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf.gz
	$tabix/tabix -p vcf $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf.gz
	$tabix/tabix -p vcf $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf.gz

	$tabix/tabix -fhB $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf.gz $onTarget_BED_169genes > $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.final.HaplotypeCaller.indel.vcf

cp $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.final.HaplotypeCaller.indel.vcf $MyOutDir/$SampleID/Target/HaplotypeCaller/$onTarget_BAM.final.HaplotypeCaller.indel.vcf

	$tabix/tabix -fhB $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf.gz $onTarget_BED_169genes > $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.final.HaplotypeCaller.snp.vcf

cp $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.final.HaplotypeCaller.snp.vcf $MyOutDir/$SampleID/Target/HaplotypeCaller/$onTarget_BAM.final.HaplotypeCaller.snp.vcf

$Java7 -$memory -jar $snpEff/SnpSift.jar tstv $MyOutDir/$SampleID/Target/HaplotypeCaller/$onTarget_BAM.final.HaplotypeCaller.snp.vcf > $MyOutDir/$SampleID/Target/HaplotypeCaller/$onTarget_BAM.final.HaplotypeCaller.snp.vcf.all.snp.TsTv.txt

	$tabix/tabix -fhB $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf.gz $onTarget_BED_169genes > $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.GVCF.vcf

cp $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.GVCF.vcf $MyOutDir/$SampleID/Target/HaplotypeCaller/$onTarget_BAM.HaplotypeCaller.snp.indel.GVCF.vcf

	$tabix/tabix -fhB $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf.gz $onTarget_BED_169genes >  $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.final.UnifiedGenotyper.indel.vcf

cp $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.final.UnifiedGenotyper.indel.vcf $MyOutDir/$SampleID/Target/UnifiedGenotyper/$onTarget_BAM.final.UnifiedGenotyper.indel.vcf

	$tabix/tabix -fhB $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf.gz $onTarget_BED_169genes > $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.final.UnifiedGenotyper.snp.vcf

cp $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.final.UnifiedGenotyper.snp.vcf $MyOutDir/$SampleID/Target/UnifiedGenotyper/$onTarget_BAM.final.UnifiedGenotyper.snp.vcf

$Java7 -$memory -jar $snpEff/SnpSift.jar tstv $MyOutDir/$SampleID/Target/UnifiedGenotyper/$onTarget_BAM.final.UnifiedGenotyper.snp.vcf > $MyOutDir/$SampleID/Target/UnifiedGenotyper/$onTarget_BAM.final.UnifiedGenotyper.snp.vcf.all.snp.TsTv.txt

############### Remove markDup, realigned BAM files
rm -rf $MyOutDir/"$SampleID"/"$SampleID".bam
rm -rf $MyOutDir/"$SampleID"/"$SampleID".bai
rm -rf $MyOutDir/"$SampleID"/"$SampleID".bam.bai
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.bai
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.bam.intervals
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bai
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam.bai
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.bam.recalibTable.grp
rm -rf $MyOutDir/"$SampleID"/"$SampleID".markDup.Realigned.recalibrated.bam.bai
#########################
	rm $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.indel.filtered.vcf.gz*
	rm $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.filtered.vcf.gz*
	rm $sftp/$RunName/HaplotypeCaller/$SampleID/$onTarget_BAM.HaplotypeCaller.snp.indel.g.vcf.gz*
	rm $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.indel.filtered.vcf.gz*
	rm $sftp/$RunName/UnifiedGenotyper/$SampleID/$onTarget_BAM.UnifiedGenotyper.snp.filtered.vcf.gz*
rm -r $sftp/$RunName/UnifiedGenotyper/H20
rm -r $sftp/$RunName/HaplotypeCaller/H20
###########################################################

# TopDir=/data/results/$Platform
TopDir=/data5/sftp/nhrisguest/incoming/upload/to_BRU

RunDate=$RunName

#bedtools=/data/Install/BEDTools-Version-2.11.2/bin/
bedtools=/data1/apps/bedtools/bedtools-2.17.0/bin/
VEP=/data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor

VEP_v=83
ASSEMBLY=GRCh37
VEP=/data1/apps/VEP/ensembl-tools-release-"$VEP_v"/scripts/variant_effect_predictor
VEP_CACHE=$HOME/.vep/
FASTA=$VEP_CACHE/homo_sapiens/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa


#LOFTEE=/data/Install/LOFTEE/loftee-master/src
LOFTEE=/data1/apps/LOFTEE/loftee-master/src
        
#MyOutDir=$TopDir/$RunDate

cd /data5/sftp/nhrisguest/incoming/upload/to_BRU/$RunDate
 		
for caller in HaplotypeCaller UnifiedGenotyper
do

	### With All samples
	MyOutDir=/data5/sftp/nhrisguest/incoming/upload/to_BRU/$RunDate

		OnTargetFile=$SampleID.markDup.Realigned.recalibrated.OnTarget.q15.bam.final
		mkdir $MyOutDir/$caller/"$SampleID"/"VEP"

	for variantType in snp indel
	do
		inVEP=$MyOutDir/$caller/"$SampleID"/$OnTargetFile.$caller.$variantType.vcf
		outVEP=$MyOutDir/$caller/"$SampleID"/VEP/$OnTargetFile.$caller.$variantType.VEP.vcf
		
		VariantCount=`grep -cv '#' $inVEP`
		
		if [ $VariantCount == 0 ]
		then
			echo -e "\nERROR : No Variants in this VCF : $inVEP\n"
			continue
		fi

		perl $VEP/variant_effect_predictor.pl \
		--dir_cache $VEP_CACHE \
		--fasta $FASTA --assembly $ASSEMBLY \
		--cache --offline --no_progress --fork $nt -i $inVEP  -o $outVEP \
		--allele_number --hgvs --check_existing --canonical --ccds --maf_exac --pubmed --protein --sift b --polyphen b \
		--fields ALLELE_NUM,Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED  \
		--vcf --force_overwrite


		# VEP QC to check if all variants are annotated. 
		# $bedtools/intersectBed -v : Only report those entries in A that have _no overlaps_ with B. - Similar to "grep -v" (an homage).

		variants_missing=`$bedtools/intersectBed -v -a $inVEP -b $outVEP | wc -l`

		if [ $variants_missing != 0 ]
		then
			echo "ERROR : These variants are missing in annotated VCF : $outVEP"
			$bedtools/intersectBed -v -a $inVEP -b $outVEP
		else
			echo -e "\nDone : VEP analysis were completed for $SampleID $caller $variantType"
		fi

		#Tablelise VEP Annotated VCFs
		out_tableize=$MyOutDir/$caller/"$SampleID"/VEP/$OnTargetFile.$caller.$variantType.VEP.txt
				
		
		python $LOFTEE/tableize_vcf.py \
		--vcf $outVEP --out $out_tableize --split_by_transcript --all_csqs --do_not_minrep --include_id \
		--info ABHet,ABHom,AC,AF,AN,DP,FS,HaplotypeScore,MLEAC,MLEAF,MQ,MQ0,QD \
		--vep_info Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED \
		--mysql

		
		ls $out_tableize
		md5sum $inVEP > $inVEP.md5
		md5sum $outVEP > $outVEP.md5
		md5sum $out_tableize > $out_tableize.md5

		echo -e "\nDone : tableize_vcf analysis were completed for $SampleID $caller $variantType"

	done
done

rm -rf $sftp/$RunName/UnifiedGenotyper/H20
rm -rf $sftp/$RunName/HaplotypeCaller/H20

date
done
