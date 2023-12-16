#!/usr/bin/bash

function get_data { fasterq-dump $1 -O $1; }

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
		for file in "${sample_data[@]}"; do if [ ! -f "$file" ]; then get_data "$file"; fi; done;
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
		for file in "${sample_data[@]}" "${control_data[@]}"; do if [ ! -f "$file" ]; then get_data "$file"; fi; done;
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

rclone copy ~/nextflow-peakcaller/output OneDrive:ChIP-seq-enrichment-analysis/Nextflow_backup
#rm -r output
