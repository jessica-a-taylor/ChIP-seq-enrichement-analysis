mv mergedBroadPeaks.bed "${sample_data[@]}"mergedBroadPeaks.bed#!/usr/bin/bash

function get_data { prefetch $1; cd $1; fastq-dump $1; rm $1.sra; cd -; }
function get_data_paired { prefetch $1; cd $1; fastq-dump $1 --split-files; rm $1.sra; cd -; }

echo "seqnamesstart, end, peak, stuff1, ., stuff2, stuff3, stuff4, counts" > headerfile

#function get_read_counts { samtools view -h $1trimmed_readsalignedReads.sam > $1trimmed_readsalignedReads.bam; samtools sort $1trimmed_readsalignedReads.bam -o $1trimmed_sorted_readsalignedReads.bam; bedtools bamtobed -i $1trimmed_sorted_readsalignedReads.bam > $1_fragments.bed; rm -rf *.sam; rm -rf *.bam; bedmap --echo --count $1__peakCallingResultsNarrow/peaks_peaks.narrowPeak $1_fragments.bed > $1_counts.bed; ssconvert $1_counts.bed $1_counts.csv; cat headerfile $1_counts.csv > $1_final_counts.csv; rm $1_counts.csv; rm $1_counts.bed; }

echo "Is a control included - yes or no?"
read control_included

echo "Sequencing type - single or paired?"
read sequencing_type

echo "Provide sample datasets"
read -a sample_data
input_sample_data=$(echo "${sample_data[*]}" | tr ' ' ,)

if [ "$control_included" = "yes" ];
then
	echo "Provide control datasets";
	read -a control_data;
	input_control_data=$(echo "${control_data[*]}" | tr ' ' ,)
fi

chmod +x bin/getFastaLength.py

if [ "$control_included" = "no" ];
then
	if [ "$sequencing_type" = "single" ];
	then
		for file in "${sample_data[@]}"; do if [ ! -f "$file" ]; then get_data "$file"; fi; done;
		nextflow run callPeaks.nf --genome-fasta TAIR10_chr_all.fas --chip-seq-fastq "$input_sample_data" --use-rmdup;
	elif [ "$sequencing_type" = "paired" ];
	then
		for file in "${sample_data[@]}"; do if [ ! -f "$file" ]; then get_data_paired "$file"; fi; done;
		nextflow run callPeaks.nf --genome-fasta TAIR10_chr_all.fas --chip-seq-fastq "$input_sample_data" --paired --use-rmdup;
	fi
elif [ "$control_included" = "yes" ];
then
	if [ "$sequencing_type" = "single" ];
	then
		for file in "${sample_data[@]}" "${control_data[@]}"; do if [ ! -f "$file" ]; then get_data "$file"; fi; done;
		nextflow run callPeaks.nf --genome-fasta TAIR10_chr_all.fas --chip-seq-fastq "$input_sample_data" --control-fastq "$input_control_data" --use-rmdup;
	elif [ "$sequencing_type" = "paired" ];
	then
		for file in "${sample_data[@]}" "${control_data[@]}"; do if [ ! -f "$file" ]; then get_data_paired "$file"; fi; done;
		nextflow run callPeaks.nf --genome-fasta TAIR10_chr_all.fas --chip-seq-fastq "$input_sample_data" --control-fastq "$input_control_data" --paired --use-rmdup;
	fi
fi

rm -rf SRR*;
cd work; rm -rf *; cd ..;
cd output;
rm -rf SRR*;
mv mergedBroadPeaks.bed "${sample_data[0]}"mergedBroadPeaks.bed;
mv mergedNarrowPeaks.bed "${sample_data[0]}"mergedNarrowPeaks.bed;
cd ..;

#for file in "${sample_data[@]}"; do get_read_counts "$file"; done;

rclone copy ~/nextflow-peakcaller/output OneDrive:ChIP-seq-enrichement-analysis/Nextflow_backup
rm -r output

