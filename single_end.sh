#!/bin/bash

#Activate the rna-seq environment
eval "$(conda shell.bash hook)"
conda activate rna-seq

# Check if 'annotation.gff' exists in the 'ref' folder
if [ -e ./ref/annotation.gff ]; then
    echo "Found annotation.gff in 'ref' folder. Converting to annotation.gtf using gffread."
    sleep 3s
    gffread ./ref/annotation.gff -T -o ./ref/annotation.gtf
elif [ -e ./ref/annotation.gtf ]; then
    echo "Found annotation.gtf in 'ref' folder. Proceeding with the existing annotation.gtf."
else
    echo "Neither annotation.gff nor annotation.gtf found in 'ref' folder."
    echo "Please provide either annotation.gff or annotation.gtf file."
fi

mkdir alignment 
mkdir FeatureCounts

for sample in $(ls -l *.fastq.gz | awk '{print $9}' | sed 's/.fastq.gz//g')
do
echo "Started running alignment using hisat2 for ${sample}"
sleep 3s
R1=${sample}.fastq.gz
echo "read 1 is ${R1}"
sleep 3s
hisat2 -x ./ref/index -U ${R1} -S ./alignment/${sample}.sam -p 60
echo "Finished running alignment using hisat2 for ${sample}"
sleep 3s
echo "Started running sam>bam conversion of alignment file using samtools for ${sample}"
sleep 3s
samtools sort ./alignment/${sample}.sam -o ./alignment/${sample}_sorted.bam
echo "Finished running sam>bam conversion of alignment file using samtools for ${sample}"
sleep 3s
echo "Started running transcript assembly using FeatureCounts for ${sample}"
sleep 3s
featureCounts -T 64 -F GTF -a ./ref/annotation.gtf -o ./FeatureCounts/${sample}_counts.txt ./alignment/${sample}_sorted.bam
echo "Finished running transcript assembly using FeatureCounts for ${sample}"
sleep 3s
done

echo "Started merging FeatureCounts output files of all samples"
sleep 3s
ls -1  FeatureCounts/*_counts.txt | parallel 'cat {} | sed '1d' | cut -f7 {} > FeatureCounts/{/.}_clean.txt'
ls -1  FeatureCounts/*_counts.txt | head -1 | xargs cut -f1 > FeatureCounts/genes.txt
paste FeatureCounts/genes.txt FeatureCounts/*_counts_clean.txt > FeatureCounts/merged_counts.txt
echo "Finished merging FeatureCounts output files of all samples"
