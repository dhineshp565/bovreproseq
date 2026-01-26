#!/usr/bin/env bash

# This script takes the sam file from minimap2,
# filters q30 and primary alignments and creates consensus of each amplicon in a sample
# generates read statistics on filtered and unfiltered read statistics
# ouput are used to detect the presence of the pathogen and the typing
# $1 = SampleName 
# $2 = Input path of sam file from minimap2, 
# $3 = primerbed file with primer coordinates for primer trimming $4 = threshold for minimum reads to be considered for consensus
# $4 = read count threshold
# $5 = samtools consensus mode "simple or bayesian"
# $6 = qscore

# generate stats prior to read filtering
samtools view -b -h $2|samtools sort > $1_unfilt.bam
samtools stats "$1_unfilt.bam" > $1_unfilt_stats.txt
samtools index "$1_unfilt.bam" > $1_unfilt.bai
samtools idxstats "$1_unfilt.bam" > $1_unfilt_idxstats.csv
#generate a  sorted bam file with primary alignments
samtools view -b -h -F 0x900 -q $6 $2|samtools sort > $1.bam	
samtools stats "$1.bam" > $1_stats.txt
#index sorted bam file and generate read counts for each amplicon
samtools index "$1.bam" > $1.bai
samtools idxstats "$1.bam" > $1_idxstats.txt
threshold=$4
awk -v var="$threshold" '{if ($3 >= var) print $1, $2, $3}' "${1}_idxstats.txt" > "${1}_mappedreads.txt"

mode=$5
#conditional for empty mapped reads.txt file
if [ $(wc -l < "$1_mappedreads.txt") -ne 0 ]
then 
#using the list of mapped amplicons from text file, bam file is split based on amplicons and consensus is generated
	while read lines
	do 
			amp=$(echo $lines|cut -f1 -d' ')
			len=$(echo $lines|cut -f2 -d' ')
			 # get upper and lower bounds of reference span

      		upper=`echo $((${len}+(${len}*20/100)))`
      		lower=`echo $((${len}-(${len}*20/100)))`

			# split bam 
			samtools view -b "$1.bam" "${amp}" > $1_${amp}.bam
			# Only reads length with + or - 10% bases is used for consenus
			samtools view -h "$1_${amp}.bam"|awk -v u=${upper} -v l=${lower} '/^@/|| length($10)>=l && length($10)<=u'|samtools sort > $1_${len}_${amp}.bam
			# generate stats for near full length reads
			samtools index "$1_${len}_${amp}.bam" > $1_${len}_${amp}.bai
			samtools idxstats "$1_${len}_${amp}.bam" > $1_${len}_${amp}_idxstats.txt
			awk '{if ($3 >= 0) print $1,$2,$3}' "$1_${len}_${amp}_idxstats.txt" > $1_${amp}_mappedreads.txt
			#trims primers from both ends of the amplicon using primer bed file
			samtools ampliconclip --both-ends -b $3 "$1_${len}_${amp}.bam" > $1_trimmed_${len}_${amp}.bam
			# generate consensus for full length reads
			if [ "$mode" == "simple" ];then
				samtools consensus -f fasta -m simple -d "$threshold" -A -q "$1_trimmed_${len}_${amp}.bam" > $1_${amp}.fasta
				# change fasta header with sample and amplicon names
				sed -i "s/>.*/>$1_${amp}_consensus/" $1_${amp}.fasta
			else
				samtools consensus -f fasta -m bayesian -A "$1_trimmed_${len}_${amp}.bam" > $1_${amp}.fasta
				# change fasta header with sample and amplicon names
				sed -i "s/>.*/>$1_${amp}_consensus/" $1_${amp}.fasta
			fi

	done < "$1_mappedreads.txt"
		# merge consensus from all amplicons
		cat $1_*.fasta > $1_cons.fasta
		cat $1_*_mappedreads.txt > $1_full_length_mappedreads.txt
		# convert multiline fasta to single line fasta
	awk '/^>/ {if (seq) print seq; print;seq =""} /^[^>]/{seq=seq$0} END {print seq}' $1_cons.fasta > $1_consensus.fasta
		
# handle empty consensus. when there are no mapped reads.add sequence header
else
		echo -e ">$1 No consensus" > $1_consensus.fasta
		echo -e "NA NA NA" > "$1_mappedreads.txt"
		echo -e "NA NA NA" > "$1_full_length_mappedreads.txt"
fi
	# insert headers to mappedreads.txt
sed -i "1i Amplicon_Name Size $1" "$1_mappedreads.txt"
sed -i "1i Amplicon_Name Size $1" "$1_full_length_mappedreads.txt"