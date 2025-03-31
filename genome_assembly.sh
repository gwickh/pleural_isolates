#!/usr/bin/env bash
#author:        :Gregory Wickham
#date:          :20250331
#version        :1.2
#date modified: :20250326
#desc           :Various tools for assembling genomes
#usage		    :bash genome_assembly.sh
#========================================================================================================

source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment
conda create -y -n genome_assembly -c bioconda -c conda-forge \
        trimmomatic shovill bowtie2 samtools bakta quast fastqc multiqc checkm2 metabat2 abricate amrfinderplus
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
    if [ -e alignment_files/${base}_alignment/${base}.bam.bai ]; then
        echo "alignment_files/${base}_alignment/${base}.bam already exists"
    else
        echo "alignment_files/${base}_alignment/${base}.sam does not exist"
        mkdir alignment_files/${base}_alignment

        bowtie2-build \
            assemblies/${base}_assembly/contigs.fa \
            alignment_files/${base}_alignment/${base}_alignment_index

        bowtie2 \
            -x alignment_files/${base}_alignment/${base}_alignment_index \
            -1 $k \
            -2 $(ls trimmed_paired/${base}_R2_*_trim.fastq.gz) \
            -S alignment_files/${base}_alignment/${base}.sam
        
        samtools sort \
            alignment_files/${base}_alignment/${base}.sam \
            -o alignment_files/${base}_alignment/${base}.bam

        samtools index ${base}.bam
    fi
done

# Calculate contig depth
for k in alignment_files/*/*.bam; do 
    jgi_summarize_bam_contig_depths \
        --outputDepth $(dirname $k)/depth.txt \
        $k
done

# Bin genomes with metabat2
mkdir genome_bins

for k in assemblies/*/*filtered.fasta; do
    dname=$(dirname $k)                             
    prefix=$(basename $dname)                       
    baseprefix=${prefix%%_assembly}              

    metabat2 \
        -a alignment_files/${baseprefix}_alignment/depth.txt \
        -i $k \
        -o genome_bins/${baseprefix}_binned \
        --minCV 0.5 \
        --minContig 1500
done

# Rename .1.fa files to .fa if .2.fa doesn't exist
for f in *.fa; do
    if [[ "$f" =~ ^(.+)\.([0-9]+)\.fa$ ]]; then 
        base="${BASH_REMATCH[1]}"     
        num="${BASH_REMATCH[2]}";   
    else     
        base="${f%.fa}";
        num=""
    fi
  # Proceed only if this is a ".1.fa" file
  if [[ "$num" == "1" ]]; then
    paired="${base}.2.fa"
    if [[ ! -f "$paired" ]]; then
    #   mv "$f" "${base}.fa"
      echo "Renamed $f → ${base}.fa"
    fi
  fi
done

# Check genome bins with checkM
mkdir checkm_reports

for k in genome_bins/*binned.fa; do
    if [ -e checkm_reports/$(basename $k)/checkm2.log ]; then
        echo "checkm_reports/$(basename $k)/checkm2.log already exists"
    else
        echo "running checkm on $k"
        checkm2 predict \
            --threads 8 \
            --input $k \
            --output-directory checkm_reports/$(basename $k)
    fi
done

(head -n 1 $(ls */*.tsv | head -n 1) 
for f in */*.tsv; do 
    sed -n '2p' "$f"; done
) > combined.tsv

# Run abricate
for k in genome_bins/*binned.fa; do
    filename=$(basename "$k")
    base=${filename%%_binned.fa}
    for j in {megares,vfdb,resfinder,card,argannot,plasmidfinder}; do
        if [ ! -f abricate/$j ]; then
            echo "making abricate/$j..."
            mkdir -p abricate/$j
        fi
        abricate \
            --db $j \
            --threads 8 \
            $k > abricate/$j/${base}_abricate_${j}.tsv
    done
done

for k in abricate/*; do
    {
        head -n 1 "$(ls $k/*.tsv | head -n 1)"
        for f in $k/*.tsv; do
            tail -n +2 "$f"
        done
    } > ${k}_combined.tsv
done

# Annotate genome with bakta
for k in assemblies/*/contigs.fa; do
    bakta \
        --db /home/ubuntu/webber_group/Gregory_Wickham/plueral_isolates/assemblies/db \
        $k \
        -o bakta_assembly/${base}_assembly.fasta \
        -t 8
done

# run refseq-masher
mkdir refseq_masher
for k in genome_bins/*.fa; do
    base=$(basename "$k" .fa)
    refseq_masher matches \
        $k \
        --output-type tab \
        -o refseq_masher/${base}_refseq_masher.tsv
done

# concatentate refseq-masher reports
(head -n 1 $(ls refseq_masher/*.tsv | head -n 1) 
for f in refseq_masher/*.tsv; do 
    sed -n '2p' "$f"; done
) > refseq_masher/combined.tsv

# clean up bakta gff3 files
for k in bakta/*/*binned.gff3; do
    dir=$(dirname "$k")
    base=$(basename "$k" .gff3)

    awk '$0 ~ /^##FASTA/ {exit} {print}' "$k" | \
    awk -F'\t' 'BEGIN { OFS="\t" }
    {
        if ($3 == "CDS" && $9 ~ /ID=/) {
            # Extract ID value
            match($9, /ID=([^;]+)/, id)
            id_val = id[1]

            # Split attributes and rebuild without ID and Name
            split($9, fields, ";")
            new_fields = ""
            for (i in fields) {
                if (fields[i] ~ /^ID=/ || fields[i] ~ /^Name=/) continue
                if (new_fields == "") {
                    new_fields = fields[i]
                } else if (fields[i] != "") {
                    new_fields = new_fields ";" fields[i]
                }
            }

            # Add new Name and optional Parent
            $9 = "Name=" id_val
            if (new_fields != "") {
                $9 = $9 ";" new_fields
            }
            if ($9 !~ /Parent=/ && id_val != "") {
                $9 = $9 ";Parent=" id_val
            }
        }
        print 
    }' > "$dir/${base}_cleaned.gff3"

    echo "✔️ Cleaned and renamed ID→Name: $dir/${base}_cleaned.gff3"
done
   

# Run AMRFinderplus
mkdir amrfinderplus

for k in bakta/*/*binned.faa; do
    amrfinder \
        -p $k \
        -n bakta/$(basename $k .faa)/$(basename $k .faa).fna \
        -g bakta/$(basename $k .faa)/$(basename $k .faa)_cleaned.gff3 \
        -o amrfinderplus/$(basename $k .faa)_amrfinderplus.tsv \
        --plus \
        --threads 8
done

(head -n 1 $(ls amrfinderplus/*amrfinderplus.tsv | head -n 1) 
for f in amrfinderplus/*amrfinderplus.tsv; do 
        tail -n +2 "$f"
done
) > amrfinderplus*amrfinderplus_combined.tsv