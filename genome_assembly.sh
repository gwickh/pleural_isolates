#!/usr/bin/env bash
#author:        :Gregory Wickham
#date:          :20250324
#version        :1.0
#date modified: :20250324
#desc           :Various tools for assembling genomes
#usage		    :bash genome_assembly.sh
#========================================================================================================

source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment
conda create -y -n genome_assembly -c bioconda -c conda-forge \
        trimmomatic shovill bowtie2 samtools bakta quast fastqc multiqc
conda activate genome_assembly

# Trim reads with trimmomatic
mkdir trimmed_paired trimmed_unpaired

for k in raw_reads/*R1_*.fastq.gz; do 
    filename=$(basename "$k")
    base=${filename%%_R1_*.fastq.gz}
    R2_file=$(ls raw_reads/${base}_R2_001.fastq.gz 2>/dev/null || ls raw_reads/${base}_R2_combined.fastq.gz 2>/dev/null)

    trimmomatic PE -threads 8 \
        "$k" "$R2_file" \
        "trimmed_paired/${base}_R1_paired_trim.fastq.gz" \
        "trimmed_unpaired/${base}_R1_unpaired_trim.fastq.gz" \
        "trimmed_paired/${base}_R2_paired_trim.fastq.gz" \
        "trimmed_unpaired/${base}_R2_unpaired_trim.fastq.gz" \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# perform read QC
mkdir fastqc_reports
fastqc trimmed_paired/* -o fastqc_reports
multiqc fastqc_reports


# Assemble genomes with spades
mkdir assemblies

for k in raw_reads/*R1_*.fastq.gz; do 
    filename=$(basename "$k")
    base=${filename%%_R1_*.fastq.gz}
    R2_file=$(ls raw_reads/${base}_R2_001.fastq.gz 2>/dev/null || ls raw_reads/${base}_R2_combined.fastq.gz 2>/dev/null)

    shovill \
        --outdir assemblies/${base}_assembly \
        --R1 "trimmed_paired/${base}_R1_paired_trim.fastq.gz"  \
        --R2 "trimmed_paired/${base}_R2_paired_trim.fastq.gz" \
        --cpus 8
done

# Align reads against assembly with bowtie2
mkdir alignment_files

for k in trimmed_paired/*_R1_*_trim.fastq.gz; do
    filename=$(basename "$k")
    base=${filename%%_R1_*_trim.fastq.gz}

    mkdir alignment_files/${base}_alignment
    bowtie2-build \
        assemblies/${base}_assembly/contigs.fa \
        alignment_files/${base}_alignment/${base}_alignment_index

    bowtie2 \
        -x alignment_files/${base}_alignment/${base}_alignment_index \
        -1 $k \
        -2 $(ls trimmed_paired/${base}_R2_*_trim.fastq.gz) \
        -S alignment_files/${base}_alignment/${base}.sam \
            | samtools sort -o ${base}.bam
    
    samtools index ${base}.bam
done

# Annotate genome with bakta
conda activate bakta

for k in assemblies/*/contigs.fa; do
    bakta \
        --db /home/ubuntu/webber_group/Gregory_Wickham/plueral_isolates/assemblies/db \
        -1 $k \
        -o bakta_assembly/${base}_assembly.fasta \
        -t 8
done
