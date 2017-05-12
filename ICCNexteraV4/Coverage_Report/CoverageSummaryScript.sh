#!/bin/bash
##Collates summary info on coverage etc for a sequencing run..

########################## Coverage Report#######################
## example to Run
## bash CoverageSummaryScript.sh 

#genomesize=3095693981	#ens60, #/data/Store/reference/human/ENSEMBL_60/fasta/dna/allchrom.validated.properties
genomesize=3137161264	#hg19, 	#from http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics

###################################### Check these variables######################
##Run details
RunDate=$1
Platform=$2
PoolName=$3

##################################################################################
logfile="$RunDate"_CoverageSummaryScript.log
exec &> >(tee -a "$logfile")

##################################################################################

# TopDir=/data/results/$Platform
TopDir=/data1/seq_data/NHCS/$Platform/results

##################################################################################

mkdir -p $TopDir/$RunDate/"Coverage_Report_v2"

#############################################################################



for RunID in $RunDate
do
	#for TargetOut in Target ProteinCodingTarget CanonTranCodingTarget
	for TargetOut in ProteinCodingTarget ProteinCodingTarget169gene CanonTranCodingTarget
	do
	
		echo "Starting Analysis for $RunID & $TargetOut"

		mkdir -p $RunID
		#Enrich_Report="Enrichment_Report"
		#TargetOut1=Target
		TargetOut1=ProteinCodingTarget	

		#dir=$TopDir/$RunID/$Enrich_Report/
		dir=$TopDir/$RunID/$PoolName/gatk_snp_indel/
		ReadsDir=$TopDir/Reads/$RunDate/
		VCFsDir=$TopDir/$RunID/$PoolName/gatk_snp_indel/
		mkdir $TopDir/$RunID/Coverage_Report_v2/"$TargetOut"
		CovDir=$TopDir/$RunID/Coverage_Report_v2/"$TargetOut"
	
		### With All samples
		SAMPLES=`ls $dir`
	
		echo -e "+++++Coverage report by $TargetOut file+++++ " 
		echo "Sample" > $RunID/Samples.tmp

		### With all samples
		ls /$dir >> $RunID/Samples.tmp

		##	Initialise temporary files
	
		echo -e "FastQC_Read1_Per_base_sequence_quality\tFastQC_Read1_Per_sequence_quality_scores" > $RunID/fastqR1.tmp
		echo -e "FastQC_Read2_Per_base_sequence_quality\tFastQC_Read2_Per_sequence_quality_scores" > $RunID/fastqR2.tmp
		echo -e "Total_HC_HetSNPs\tHetSNPs_AB<0.40or>0.60\t%HetSNPs_AB<0.40or>0.60" > $RunID/HetSNPs.tmp
		echo -e "TotalReads\tMappedReads\t%Mapped\tMappedReads_q15\t%Mapped_q15\tReadsOnDesign_q15\t%OnDesign_q15\tReadsOnTarget_q15\t%OnTarget\tUniqReadsOnTarget_q15\t%UniqueReadsOnTarget_q15" > $RunID/Reads.tmp

		echo -e "MeanFwdReadLength\tMeanRevReadLength\tReadsAlignedInPairs\t%ReadsAlignedInPairs\tStrandBalance" > $RunID/Reads2.tmp
		echo "EnrichmentFactor" > $RunID/EF.tmp
		echo "ExonsMissed" > $RunID/ExonsMissed.tmp
		echo -e "Bases=0x\tBases>=1x\tBases>=5x\tBases>=10x\tBases>=20x\tBases>=30x\tBases>=40x\tBases>=50x" > $RunID/Coverage.tmp
		echo -e "TargetSize(bp)\tMeanCov\tMaxCov\tMedianCov" > $RunID/CoverageStats.tmp
		echo -e "Callable\t%Callable\tNo_Coverage\tLow_Coverage\tExcess_Coverage\tPoor_Quality\tTotal" > $RunID/Callable.tmp
		echo -e "MeanCov\tMedianCov\tDiff(mean-med)\tMaxCov\tTarget_size\tEvenness" > $RunID/Evenness.tmp
		echo -e "SNPs_Ts/Tv_UnifiedGenotyper" > $RunID/TsTvUG.tmp
		echo -e "SNPs_Ts/Tv_HaplotyeCaller" > $RunID/TsTvHC.tmp
		#####
		echo -e "tabulating read data & calculating enrichment factor"
		
		for sample in $SAMPLES
		do
			echo -n "$sample "
			##	Counts total reads / mapped reads (Total, Mapped (q8), % Mapped, OnDESIGN, OnTARGET) 
			TotR=`cat $dir$sample/*.TotalReads.flagstat| awk 'NR==1 {print $1}'`
			TotR_mapped=`cat $dir$sample/*.mapped.flagstat | awk 'NR==5 {print $1}'`  ## old samtools NR==4 into new NR==5
			TotR_mapped_pc=`echo $TotR $TotR_mapped | awk '{printf "%.2f\n", 100*$2/$1}'`
			TotR_mapped_q15=`cat $dir$sample/*.Realigned.recalibrated.bam.flagstat  | awk 'NR==5 {print $1}'`  ##NR==4 into NR==5
			TotR_mapped_q15_pc=`echo $TotR_mapped $TotR_mapped_q15 | awk '{printf "%.2f\n", 100*$2/$1}'`
			TotR_ondesign_pc=`echo $TotR_mapped_q15 $TotR_ondesign | awk '{printf "%.2f\n", 100*$2/$1}'`
			TotR_ontarget=`cat $dir$sample/$TargetOut/*.flagstat | awk 'NR==1 {print $1}'`
			TotR_ontarget_pc=`echo $TotR_mapped_q15 $TotR_ontarget | awk '{printf "%.2f\n", 100*$2/$1}'`
			DupR_ontarget=`cat $dir$sample/$TargetOut/*.bam.flagstat | awk 'NR==4 {print $1}'`  ##NR==3 into NR==2
			DupR_ontarget_pc=`echo $TotR_ontarget $DupR_ontarget | awk '{printf "%.2f\n", 100*$2/$1}'`
			UniqR_ontarget=`echo $TotR_ontarget $DupR_ontarget | awk '{print $1-$2}'`
			UniqR_ontarget_pc=`echo $TotR_ontarget $UniqR_ontarget | awk '{printf "%.2f\n", 100*$2/$1}'`


			echo -e "$TotR\t$TotR_mapped\t$TotR_mapped_pc\t$TotR_mapped_q15\t$TotR_mapped_q15_pc\t$TotR_ondesign\t$TotR_ondesign_pc\t$TotR_ontarget\t$TotR_ontarget_pc\t$UniqR_ontarget\t$UniqR_ontarget_pc" >> $RunID/Reads.tmp

			##	Read lengths & Strand balance (Only considers reads mapping to target)
			ReadL_F=`cat $dir$sample/$TargetOut/*.AlignSumMet.txt | awk '$1 !~ /^#/ && NF !=0' | awk '{FS="\t"; if (NR==1) for (i = 1; i <= NF; i++) if ($i=="MEAN_READ_LENGTH") n=i ; print $n}' | awk 'NR==2' | awk '{printf "%.2f\n",$1}'`
			ReadL_R=`cat $dir$sample/$TargetOut/*.AlignSumMet.txt | awk '$1 !~ /^#/ && NF !=0' | awk '{FS="\t"; if (NR==1) for (i = 1; i <= NF; i++) if ($i=="MEAN_READ_LENGTH") n=i ; print $n}' | awk 'NR==3' | awk '{printf "%.2f\n",$1}'`
			Reads_aligned_in_pairs=`cat $dir$sample/$TargetOut/*.AlignSumMet.txt | awk '$1 !~ /^#/ && NF !=0' | awk '{FS="\t"; if (NR==1) for (i = 1; i <= NF; i++) if ($i=="READS_ALIGNED_IN_PAIRS") n=i ; print $n}' | awk 'NR==4'`
			Reads_aligned_in_pairs_pc=`cat $dir$sample/$TargetOut/*.AlignSumMet.txt | awk '$1 !~ /^#/ && NF !=0' | awk '{FS="\t"; if (NR==1) for (i = 1; i <= NF; i++) if ($i=="PCT_READS_ALIGNED_IN_PAIRS") n=i ; print $n}' | awk 'NR==4 {printf "%.2f\n", 100*$1}'`
			StrandBalance=`cat $dir$sample/$TargetOut/*.AlignSumMet.txt | awk '$1 !~ /^#/ && NF !=0' | awk '{FS="\t"; if (NR==1) for (i = 1; i <= NF; i++) if ($i=="STRAND_BALANCE") n=i ; print $n}' | awk 'NR==4' | awk '{printf "%.3f\n",$1}'`
			echo -e "$ReadL_F\t$ReadL_R\t$Reads_aligned_in_pairs\t$Reads_aligned_in_pairs_pc\t$StrandBalance" >>$RunID/Reads2.tmp
			##	Enrichment factor
			length=`cat $dir$sample/$TargetOut/*.CoveragePerBase | wc -l`
			echo $TotR_ontarget $TotR_mapped_q15 $length $genomesize | awk '{print ($1/$2)/($3/$4)}' >> $RunID/EF.tmp
		done

		#####
		##	Counts missing exons (in TARGET interval)
		echo -e "\n Tabulating missing exons"
		for sample in $SAMPLES
		do
			echo -n "$sample "
			ExonsMissed=`awk '$5==0' $dir$sample/$TargetOut/*.CoverageStats | wc -l`
			echo $ExonsMissed >> $RunID/ExonsMissed.tmp
		done

		#####
		##	Coverage thresholds (No Coverage, 0x, >=1x, >=5x,>=10x)(relative to TARGET interval)
		echo -e "\n Tabulating coverage"
		for sample in $SAMPLES
		do
			echo -n "$sample "
		
			# (No Coverage, 0x, >=1x, >=5x,>=10x,>=20x,>=30x) by percentage
			#awk '$1 ~/all/ {if ($2==0) a=a+$3;if($2>=1) b=b+$3; if($2>=5) c=c+$3; if ($2>=10) d=d+$3; if ($2>=20) e=e+$3; if ($2>=30) f=f+$3}END{OFS="\t";print 100*a/$4,100*b/$4,100*c/$4,100*d/$4,100*e/$4,100*f/$4}' $dir$sample/$TargetOut/*.CoverageHist >> $RunID/Coverage.tmp

			## (No Coverage, 0x, >=1x, >=5x,>=10x, >=20x, >=30x, >=40x, >=50x) by percentage
			awk '$1 ~/all/ {if ($2==0) a=a+$3;if($2>=1) b=b+$3; if($2>=5) c=c+$3; if ($2>=10) d=d+$3; if ($2>=20) e=e+$3; if ($2>=30) f=f+$3;  if ($2>=40) g=g+$3; if ($2>=50) h=h+$3}END{OFS="\t";printf "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", 100*a/$4,100*b/$4,100*c/$4,100*d/$4,100*e/$4,100*f/$4,100*g/$4,100*h/$4}' $dir$sample/$TargetOut/*.CoverageHist >> $RunID/Coverage.tmp

			## (No Coverage, 0x, >=1x, >=5x,>=10x, >=20x, >=30x, >=40x, >=50x) by actual numbers
			##awk '$1 ~/all/ {if ($2==0) a=a+$3;if($2>=1) b=b+$3; if($2>=5) c=c+$3; if ($2>=10) d=d+$3; if ($2>=20) e=e+$3; if ($2>=30) f=f+$3;  if ($2>=40) g=g+$3; if ($2>=50) h=h+$3}END{OFS="\t";print a,b,c,d,e,f,g,h}' $dir$sample/$TargetOut/*.CoverageHist >> $RunID/Coverage.tmp


		done

		#####
		##	Coverage stats (relative to TARGET) 
		echo -e "\ncalculating coverage stats"
		for sample in $SAMPLES
		do
			echo -n "$sample "
			stats=`awk '{t=t+$6; if ($6>n) n=$6}END{OFS="\t";print NR,t/NR,n}' $dir$sample/$TargetOut/*.CoveragePerBase` #calculates length, mean & max coverage in 1 line
			cat $dir$sample/$TargetOut/*.CoveragePerBase | cut -f6 | sort -n > $RunID/covtemp.tmp
			length=`cat $RunID/covtemp.tmp | wc -l`
			median=`awk "BEGIN{n=$length}"'NR==n/2' $RunID/covtemp.tmp`
			echo -e "$stats\t$median" >> $RunID/CoverageStats.tmp
		done

		#####
		##	Callable Bases (relative to TARGET interval)
		echo -e "\ntabulating callable bases"
		for sample in $SAMPLES
		do
			echo -n "$sample "
			cat $dir$sample/$TargetOut/*.callable.summary | awk '{if ($1 ~ /CALLABLE/) a=$2 ; if ( $1 ~ /NO_C/) b=$2 ; if ($1 ~/LOW_C/) c=$2 ;if ($1 ~ /EXC/) d=$2 ; if ($1~ /POOR/) e=$2 ; }END{OFS="\t"; printf "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", a,100*a/(a+b+c+d+e),b,c,d,e,a+b+c+d+e}' >> $RunID/Callable.tmp
		done

		## Calculates Evenness (and med/mean/max) in R
		echo -e "\ncalculating evenness"
		for sample in $SAMPLES
		do
			CoveragePerBase=`ls $dir$sample/$TargetOut/ | grep q15.bam.CoveragePerBase` 
		echo -n "$sample "
		echo "mydata <- read.delim(\"$dir$sample/$TargetOut/$CoveragePerBase\", header=F)
				mymean <- mean(mydata\$V6)
				mymedian <- median(mydata\$V6)
				mymax <- max(mydata\$V6)
				mydiff <- mymean - mymedian

				myroundmean <- trunc(mymean)
				mytarget <- length(mydata\$V6)

				n <- rep(0,mymax)
				for (i in 1:mymax) {
					n[i] <- length(mydata\$V6[mydata\$V6>=i])
				}
				a <- sum(n[1:myroundmean])
				b <- sum(n[(myroundmean+1):mymax])
				myevenness <- 100*a/(a+b) #Evenness score as defined Mokry et al

				# %bases covered with depth 1 < x < meandepth
				c <- 100 * length(mydata\$V6[mydata\$V6!=0 & mydata\$V6<mymean]) / length(mydata$V6)

				# plot(100*n/(mytarget),type=\"l\",ylim=c(0,100))
				# abline(v=mymean)
				# abline(v=mymedian)
				write(c(mymean,mymedian,mydiff,mymax,mytarget,myevenness),ncolumns=6,file=\"$RunID/covstats.tmp\",sep=\"\\t\")
			" > $RunID/script.r

		/data1/apps/R/R-3.0.2/bin/R CMD BATCH $RunID/script.r
		cat $RunID/covstats.tmp | awk 'BEGIN {FS="\t";} {printf "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", $1,$2,$3,$4,$5,$6}'  >> $RunID/Evenness.tmp

		done

		echo -e "\nCalculating Ts/Tv ratio"
		### Ts/Tv ratio

		for sample in $SAMPLES
		do  
			echo -n "$sample "
			tstvUG=`cat $dir$sample/Target/UnifiedGenotyper/*.final.UnifiedGenotyper.snp.vcf.all.snp.TsTv.txt | awk 'NR==6' | awk 'BEGIN {FS=",";} {printf "%.2f\n", $2}'`
			echo -e "$tstvUG" >> $RunID/TsTvUG.tmp
			
			if [ "$RunID" = "$RunDate" ]
			then
				tstvHC=`cat $dir$sample/Target/HaplotypeCaller/*.final.HaplotypeCaller.snp.vcf.all.snp.TsTv.txt | awk 'NR==6' | awk 'BEGIN {FS=",";} {printf "%.2f\n",  $2}'`
				echo -e "$tstvHC" >> $RunID/TsTvHC.tmp
			fi 
		done

		echo -e "\nGetting FastQC results"
		### FastQC values Read1
		for sample in $SAMPLES
		do
			echo -n "$sample "
			read1=`ls $ReadsDir/$sample/ | grep Trimed_1_fastqc | head -1`	
			R1_fastqPerBaseQua=`cat $ReadsDir$sample/$read1/summary.txt | awk 'NR==2 {print $1}'`
			R1_fastqPerSeqQuaScore=`cat $ReadsDir$sample/$read1/summary.txt | awk 'NR==3 {print $1}'`
			echo -e "$R1_fastqPerBaseQua\t$R1_fastqPerSeqQuaScore" >> $RunID/fastqR1.tmp
		done
		### FastQC values Read2
		for sample in $SAMPLES
		do
			read2=`ls $ReadsDir/$sample/ | grep Trimed_2_fastqc | head -1`
			R2_fastqPerBaseQua=`cat $ReadsDir$sample/$read2/summary.txt | awk 'NR==2 {print $1}'`
			R2_fastqPerSeqQuaScore=`cat $ReadsDir$sample/$read2/summary.txt | awk 'NR==3 {print $1}'`
			echo -e "$R2_fastqPerBaseQua\t$R2_fastqPerSeqQuaScore" >> $RunID/fastqR2.tmp
		done
		
		
		#HaplotypeCaller Heterozygous SNPs ratio
		echo -e "\nCalculating HaplotypeCaller Heterozygous SNPs ratio"
		for sample in $SAMPLES
		do  
			echo -n "$sample "
			# Total Heterozygous SNPs
			TotHet=`cat $VCFsDir/$sample/Target/HaplotypeCaller/*.final.HaplotypeCaller.snp.vcf |  grep ABHet= | grep PASS | cut -f8 | awk 'BEGIN {FS=";"} {print $1}' | awk 'BEGIN {FS="="} {print $2}' | wc -l` 
			# Inside >=0.40 && <= 0.60
			InSideHet=`cat $VCFsDir/$sample/Target/HaplotypeCaller/*.final.HaplotypeCaller.snp.vcf | grep ABHet= | grep PASS | cut -f8 | awk 'BEGIN {FS=";"} {print $1}' | awk 'BEGIN {FS="="} {print $2}' | awk '{ if ($1 >= 0.40 && $1 <= 0.60) print $1 }' | wc -l`
			# Outside <0.40 || > 0.60
			OutSideHet=`cat $VCFsDir/$sample/Target/HaplotypeCaller/*.final.HaplotypeCaller.snp.vcf | grep ABHet= | grep PASS | cut -f8 | awk 'BEGIN {FS=";"} {print $1}' | awk 'BEGIN {FS="="} {print $2}' | awk '{ if ($1 < 0.40 || $1 > 0.60) print $1 }' | wc -l`
			
			InHetRatio=`echo $TotHet $InSideHet | awk '{printf "%.2f\n", 100*$2/$1}'`
			OutHetRatio=`echo $TotHet $OutSideHet | awk '{printf "%.2f\n", 100*$2/$1}'`
			echo -e "$TotHet\t$OutSideHet\t$OutHetRatio" >> $RunID/HetSNPs.tmp
		done

		echo -e "\nAssembles output file"
		##	Assembles output file
		paste $RunID/Samples.tmp $RunID/Reads.tmp $RunID/Reads2.tmp $RunID/EF.tmp $RunID/ExonsMissed.tmp $RunID/Coverage.tmp $RunID/CoverageStats.tmp $RunID/Callable.tmp $RunID/Evenness.tmp $RunID/TsTvUG.tmp $RunID/TsTvHC.tmp $RunID/fastqR1.tmp $RunID/fastqR2.tmp $RunID/HetSNPs.tmp > $RunID/tmp10.tmp

		## column 25,26,27 : >=30x, >=40x, >=50x
		cat $RunID/tmp10.tmp | awk 'BEGIN {FS=OFS="\t";} {print $1,$2,$3,$4,$5,$6,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$27,$28,$32,$33,$34,$35,$36,$37,$39,$40,$41,$42,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54}' > SummaryOutput_"$TargetOut"_"$RunID".txt

		## Remove tmp files
		rm $RunID/*.tmp
		rm $RunID/*.r
		rm script.r.Rout 
		rm -rf $RunID
		mv SummaryOutput_"$TargetOut"_"$RunID".txt $CovDir/
		
	done
done

echo -e "Getting FastQC summary"
## Getting FastQC summary 
bash /data1/seq_data/Store/Scripts/Research/FastQC_Report_v2.sh $RunDate $Platform $PoolName
echo -e "\n ++++++++DONE++++++++ \n"		
echo "End of Analysis for $RunID & $TargetOut"
echo -e "End Time : `date`"
exit
