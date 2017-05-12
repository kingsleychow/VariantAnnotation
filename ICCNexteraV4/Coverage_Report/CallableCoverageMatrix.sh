#!/bin/bash

## Creates matrix: % bases with callable coverage by sample and exon
## Example to Run
## bash CallableCoverageMatrix.sh 

logfile=CallableCoverageMatrix.log
exec &> >(tee -a "$logfile")

###################################### Check these variables######################
##Run details
RunDate=$1
Platform=$2
PoolName=$3

##################################################################################

# TopDir=/data/results/$Platform
TopDir=/data1/seq_data/NHCS/$Platform/results

##################################################################################

mkdir -p $TopDir/$RunDate/"Coverage_Report_v2"

echo -e "\n Start time : `date`"

# Specify the Coverage report from ProteinCodingTarget, CanonTranCodingTarget

for RunID in $RunDate
do

	echo "Starting Analysis for $RunID"

	for TargetOut in ProteinCodingTarget ProteinCodingTarget_20x ProteinCodingTarget_30x CanonTranCodingTarget CanonTranCodingTarget_20x CanonTranCodingTarget_30x ProteinCodingTarget169gene ProteinCodingTarget169gene_20x ProteinCodingTarget169gene_30x
	do 
		echo -e "\n++++CallableLoci report by : $TargetOut " 

		#dir=$TopDir/$RunID/$Enrich_Report/
		dir=$TopDir/$RunID/$PoolName/gatk_snp_indel/
		rundir=$TopDir/$RunID/$PoolName/gatk_snp_indel/
		mkdir $TopDir/$RunID/Coverage_Report_v2/"$TargetOut"
		CovDir=$TopDir/$RunID/Coverage_Report_v2/"$TargetOut"


		SAMPLES=`ls /$rundir`
		FIRSTSAMPLE=`ls /$rundir | head -1`


		cp $dir/$FIRSTSAMPLE/$TargetOut/*bam.CoverageStats stat.tmp

		sed -i 's/_/|/g' stat.tmp
		cat stat.tmp| cut -f4|  awk 'BEGIN {FS="|"} {print $1}' > stat2.tmp
		paste <(cut -f1-3 stat.tmp) <(cut -f1 stat2.tmp) <(cut -f5-8 stat.tmp) > stat3.tmp
		mv stat3.tmp temp.tmp


		AMPLICONS=`cat temp.tmp | cut -f4 | sort -u`
		##	List of intervals:
		echo "Defining Samples"
		echo -e "Chr\tStart\tEnd\tGeneID\tExonSize" > amplicons.tmp

		for amplicon in $AMPLICONS
			do
			cat temp.tmp | grep -w $amplicon | awk '{OFS="\t"; print $1,$2,$3,$4,($3-$2)}' >> amplicons.tmp
		done

		##  Find zero intervals by amplicon for each sample in turn
		for sample in $SAMPLES
		do
			echo "Processing Sample $sample"
			echo "$sample" > $sample.tmp
			cp $dir/$sample/$TargetOut/*callable.byTarget temp.tmp
		
			sed -i 's/_/|/g' temp.tmp
			cat temp.tmp | cut -f4|  awk 'BEGIN {FS="|"} {print $1}' > tar.tmp
			paste <(cut -f1-3 temp.tmp) <(cut -f1 tar.tmp) <(cut -f5-7 temp.tmp) > tar2.tmp
			mv tar2.tmp temp.tmp

			for amplicon in $AMPLICONS
				do
				cat temp.tmp | grep -w $amplicon | awk '{print 100*$5/($5+$6+$7)}' >> $sample.tmp
			done
		done


		##	Build output file by pasting each results for each sample together
		cp amplicons.tmp temp1.tmp
		echo "writing file"
		for sample in $SAMPLES
		do
			paste temp1.tmp $sample.tmp > temp2.tmp
			cp temp2.tmp temp1.tmp
		done

		mv temp1.tmp CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt 


		### Remove duplicate lines
		FileCall=CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt
		cat $FileCall | awk 'NR==1' > head.tmp
		num=`sed '1d' $FileCall |wc -l`
		echo -e "Number of Lines in PerCallableByTarget File : $num"
		seq $num > testNum.tmp
		sed '1d' $FileCall |awk 'BEGIN {FS=OFS="\t";} {print $1"_"$2"_"$3"_"$4,$0}' > PerCall1.tmp
		paste <(cut -f1 testNum.tmp) <(cut -f1- PerCall1.tmp ) > PerCall2.tmp
		awk '!array[$2]++' PerCall2.tmp > PerCall3.tmp
		sort -k1n PerCall3.tmp | cut -f3- > PerCall4.tmp
		cat head.tmp PerCall4.tmp > finalout.txt
		rm -rf *.tmp
		mv finalout.txt $FileCall
		NumLine=`wc -l $FileCall`
		echo -e "Actual number of Target including header : $NumLine"
		echo "Writing file : CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt"	

		Genes=`awk '{print $4}' CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt |sed '1d'|  sort -u`
		echo -e "Gene\tMaxCallable\tMinCallable\tAvgCallable\tStdevCallable" > CallableByRun_"$TargetOut"_"$RunID"_"$PoolName".txt

		head1=`head -n 1 CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt| paste <(cut -f4,6-)`
		#echo -e "Header1 : $head1"
		echo  -e "$head1" > CallableBySample_"$TargetOut"_"$RunID"_"$PoolName".txt

		for Gene in $Genes
		do
			## Remove Negative Control
			awk 'BEGIN {FS=OFS="\t";}
				NR==FNR{
					A[$1]
					next
					}
				{
				s=$1
			for(i=2; i<=NF; i++){
				if(FNR==1) for(j in A) if($i~j) D[i]
				if( ! (i in D) ) s=s OFS $i
				}
					print s
			}' /data1/seq_data/Store/Scripts/NegControl.File.with.patterns.txt CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt > miscPerCall

		cat miscPerCall | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc2
		cat miscPerCall | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i*$5}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc3
		paste <(cut -f1 misc2) <(cut -f2- misc3) > misc4
		
		cat CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc2a
		cat CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i*$5}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc3a
			paste <(cut -f1 misc2a) <(cut -f2- misc3a) > misc4a	

		Calla=`cat misc4a | awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i/$1}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=1`
		Call=`cat misc4 | awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i/$1}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=1`
		
		## Callable by Run -- maximum,minimum,avg,stdev 
		MaxCallable=`echo $Call | awk 'BEGIN {FS=OFS=" ";} {for (i=1;i<=NF;i++) { arr[NR,i]=$i;if(big <= NF)   big=NF; } } END { for(i=1;i<=big;i++) {for(j=1;j<=NR;j++) { printf("%s\t",arr[j,i]); } printf("\n"); } }' | sort -k1n| tail -n 1`
		
		MinCallable=`echo $Call | awk 'BEGIN {FS=OFS=" ";} {for (i=1;i<=NF;i++) { arr[NR,i]=$i;if(big <= NF)   big=NF; } } END { for(i=1;i<=big;i++) {for(j=1;j<=NR;j++) { printf("%s\t",arr[j,i]); } printf("\n"); } }' | sort -k1n| head -n 1`
		
		AvgCallable=`echo $Call | awk 'BEGIN {FS=OFS=" ";} { A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; print A/=NF }'`
		StdevCallable=`echo $Call | awk 'BEGIN {FS=OFS=" ";} { A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ; for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }'`
		
		rm -rf misc*
		rm -rf *.tmp

		echo -e "$Gene\t$Calla" >> CallableBySample_"$TargetOut"_"$RunID"_"$PoolName".txt
		echo -e "$Gene\t$MaxCallable$MinCallable$AvgCallable\t$StdevCallable" >> CallableByRun_"$TargetOut"_"$RunID"_"$PoolName".txt

		done

		mv CallableBySample_"$TargetOut"_"$RunID"_"$PoolName".txt $CovDir/
		mv CallableByRun_"$TargetOut"_"$RunID"_"$PoolName".txt $CovDir/
		mv CallableByExon_"$TargetOut"_"$RunID"_"$PoolName".txt $CovDir/
		echo -e "Callable report Completed : $TargetOut"
	done



	##############################################Mean Coverage by Exon#################################################################
	#####TargetOut=ProteinCodingTarget CanonTranCodingTarget CoreGeneCanonTranCodingTarget

	echo -e "\n++++ Mean Coverage Calculation by ProteinCodingTarget,  CanonTranCodingTarget  Exons +++++"
	for TargetOut in ProteinCodingTarget CanonTranCodingTarget 
	do 

		echo -e "Mean Coverage report by : $TargetOut "
		##dir=$TopDir/$RunID/Enrichment_Report/
		dir=$TopDir/$RunID/$PoolName/gatk_snp_indel/
		rundir=$TopDir/$RunID/$PoolName/gatk_snp_indel/
		mkdir $TopDir/$RunID/Coverage_Report_v2/"$TargetOut"
		CovDir=$TopDir/$RunID/Coverage_Report_v2/"$TargetOut"

		SAMPLES=`ls /$rundir`
		FIRSTSAMPLE=`ls /$rundir | head -1`


		cp $dir/$FIRSTSAMPLE/$TargetOut/*bam.CoverageStats stat.tmp

		sed -i 's/_/|/g' stat.tmp
		cat stat.tmp| cut -f4|  awk 'BEGIN {FS="|"} {print $1}' > stat2.tmp
		paste <(cut -f1-3 stat.tmp) <(cut -f1 stat2.tmp) <(cut -f5-8 stat.tmp) > stat3.tmp
		mv stat3.tmp temp.tmp


		AMPLICONS=`cat temp.tmp | cut -f4 | sort -u`
		##      	List of intervals:
		echo "Defining Samples"
		echo -e "Chr\tStart\tEnd\tGeneID\tExonSize" > amplicons.tmp

		for amplicon in $AMPLICONS
			do
				cat temp.tmp | grep -w $amplicon | awk '{OFS="\t"; print $1,$2,$3,$4,($3-$2)}' >> amplicons.tmp
			done

		##  Find zero intervals by amplicon for each sample in turn
		for sample in $SAMPLES
		do
			echo "Processing Sample $sample"
				echo "$sample" > $sample.tmp
				cp $dir$sample/$TargetOut/*sample_interval_summary.MeanCvg.bed temp.tmp
			
				sed -i 's/_/|/g' temp.tmp
				cat temp.tmp | cut -f4|  awk 'BEGIN {FS="|"} {print $1}' > tar.tmp
				paste <(cut -f1-3 temp.tmp) <(cut -f1 tar.tmp) <(cut -f5-7 temp.tmp) > tar2.tmp
				mv tar2.tmp temp.tmp

				for amplicon in $AMPLICONS
				do
				cat temp.tmp | grep -w $amplicon | awk '{print $5}' >> $sample.tmp
				done
		done


		##      Build output file by pasting each results for each sample together
		cp amplicons.tmp temp1.tmp
		echo "writing file"
		for sample in $SAMPLES
			do
			paste temp1.tmp $sample.tmp > temp2.tmp
				cp temp2.tmp temp1.tmp
			done

		mv temp1.tmp AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt 

		### Remove duplicate lines
			FileCall=AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt
			cat $FileCall | awk 'NR==1' > head.tmp
			num=`sed '1d' $FileCall |wc -l`
			echo -e "Number of Lines in PerCallableByTarget File : $num"
			seq $num > testNum.tmp
			sed '1d' $FileCall |awk 'BEGIN {FS=OFS="\t";} {print $1"_"$2"_"$3"_"$4,$0}' > PerCall1.tmp
			paste <(cut -f1 testNum.tmp) <(cut -f1- PerCall1.tmp ) > PerCall2.tmp
			awk '!array[$2]++' PerCall2.tmp > PerCall3.tmp
			sort -k1n PerCall3.tmp | cut -f3- > PerCall4.tmp
			cat head.tmp PerCall4.tmp > finalout.txt
			rm -rf *.tmp
			mv finalout.txt $FileCall
			NumLine=`wc -l $FileCall`
			echo -e "Actual Number of number of Target : $NumLine"

		Genes=`awk '{print $4}' AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt |sed '1d'|  sort -u`
		echo -e "Gene\tMaxCoverage\tMinCoverage\tAvgCoverage\tStdevCoverage" > AvgCoverageByRun_"$TargetOut"_"$RunID"_"$PoolName".txt

		head1=`head -n 1 AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt| paste <(cut -f4,6-)`
		#echo -e "Header1 : $head1"
		echo  -e "$head1" > AvgCoverageBySample_"$TargetOut"_"$RunID"_"$PoolName".txt

		for Gene in $Genes
		do
				## Remove Negative Control
				awk 'BEGIN {FS=OFS="\t";}
				NR==FNR{
						A[$1]
						next
						}
				{
						s=$1
						for(i=2; i<=NF; i++){
						if(FNR==1) for(j in A) if($i~j) D[i]
						if( ! (i in D) ) s=s OFS $i
				}
				print s
			}' /data1/seq_data/Store/Scripts/NegControl.File.with.patterns.txt AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt > miscPerCall

				cat miscPerCall | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc2
				cat miscPerCall | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i*$5}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc3
				paste <(cut -f1 misc2) <(cut -f2- misc3) > misc4

				cat AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc2a
				cat AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt | grep -w $Gene |  awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i*$5}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=4  > misc3a
				paste <(cut -f1 misc2a) <(cut -f2- misc3a) > misc4a

				Calla=`cat misc4a | awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i/$1}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=1`
				Call=`cat misc4 | awk 'BEGIN {FS=OFS="\t";} {for(i=c;++i<=NF;)s[i]+=$i/$1}END{for(i=c;++i in s;)printf("%.2f%c",s[i],((i+1) in s)?FS:RS)}' c=1`

				## Coverage by Run -- max,min,avg,stdev
				MaxCoverage=`echo $Call | awk 'BEGIN {FS=OFS=" ";} {for (i=1;i<=NF;i++) { arr[NR,i]=$i;if(big <= NF)   big=NF; } } END { for(i=1;i<=big;i++) {for(j=1;j<=NR;j++) { printf("%s\t",arr[j,i]); } printf("\n"); } }' | sort -k1n| tail -n 1`
			MinCoverage=`echo $Call | awk 'BEGIN {FS=OFS=" ";} {for (i=1;i<=NF;i++) { arr[NR,i]=$i;if(big <= NF)   big=NF; } } END { for(i=1;i<=big;i++) {for(j=1;j<=NR;j++) { printf("%s\t",arr[j,i]); } printf("\n"); } }' | sort -k1n| head -n 1`
			AvgCoverage=`echo $Call | awk 'BEGIN {FS=OFS=" ";} { A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; print A/=NF }'`
			StdevCoverage=`echo $Call | awk 'BEGIN {FS=OFS=" ";} { A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ; for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }'`

				rm -rf misc*
				rm -rf *.tmp

				echo -e "$Gene\t$Calla" >> AvgCoverageBySample_"$TargetOut"_"$RunID"_"$PoolName".txt
				echo -e "$Gene\t$MaxCoverage$MinCoverage$AvgCoverage\t$StdevCoverage" >> AvgCoverageByRun_"$TargetOut"_"$RunID"_"$PoolName".txt

		done

		mv AvgCoverageBySample_"$TargetOut"_"$RunID"_"$PoolName".txt $CovDir/
		mv AvgCoverageByRun_"$TargetOut"_"$RunID"_"$PoolName".txt $CovDir/
		mv AvgCoverageByExon_"$TargetOut"_"$RunID"_"$PoolName".txt $CovDir/
		echo -e "Mean Coverage report Completed : $TargetOut"
	done
	echo "End of Analysis for $RunID"

done
echo -e "+++++ALL DONE+++++"
echo -e "End time : `date`"
exit

