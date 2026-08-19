#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="ref/MN908947.3.fasta"
SAMPLES_FILE="acc_list.txt"
METADATA="metadata.csv"

RAW_DIR="raw"
TRIM_DIR="trimmed"
ALIGN_DIR="align"
CONSENSUS_DIR="consensus"
LINEAGE_DIR="lineage"
QC_DIR="qc"

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$CONSENSUS_DIR" "$LINEAGE_DIR" "$QC_DIR"

if [[ ! -f "${REF}.bwt" ]]; then
    echo "Building BWA index for REF..."
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "Indexing reference FASTA..."
    samtools faidx "$REF"
fi

echo "Updating Nextclade SARS-CoV-2 dataset..."
mkdir -p "ref/nextclade_sars-cov-2"
nextclade dataset get --name sars-cov-2 --output-dir "ref/nextclade_sars-cov-2" 

get_primer_scheme () {
    local acc="$1"
    awk -F',' -v acc="$acc" '$1 == acc { print $3 }' "$METADATA"
}

while read -r ACC || [[ -n "$ACC" ]]; do
    ACC=$(echo "$ACC" | tr -d '\r\n ')
    [[ -z "$ACC" || "$ACC" =~ ^# ]] && continue

    # SKIP IF ALREADY COMPLETE
    if [[ -f "${LINEAGE_DIR}/${ACC}.pangolin.csv" ]]; then
        echo "[$ACC] Lineage assignment already complete. Skipping..."
        continue
    fi

    echo "=================================================================="

    # Download (Prefetch + fasterq-dump)
    R1="${RAW_DIR}/${ACC}_1.fastq.gz"
    R2="${RAW_DIR}/${ACC}_2.fastq.gz"

    if [[ ! -f "$R1" && ! -f "$R2" && ! -f "${RAW_DIR}/${ACC}.fastq.gz" ]]; then
        echo "[$ACC] Raw FASTQ missing. Downloading from SRA..."
        
        # Prefetch container first to prevent network drops
        prefetch "$ACC" --output-directory "$RAW_DIR"
        
        # Extract locally
        if [[ -d "${RAW_DIR}/${ACC}" ]]; then
            fasterq-dump "${RAW_DIR}/${ACC}" --outdir "$RAW_DIR" --split-3 -e "$THREADS" --temp "$RAW_DIR"
            rm -rf "${RAW_DIR}/${ACC}"
        else
            fasterq-dump "$ACC" --outdir "$RAW_DIR" --split-3 -e "$THREADS" --temp "$RAW_DIR"
        fi

        # Fast parallel compression
        if [[ -f "${RAW_DIR}/${ACC}_1.fastq" ]]; then
            pigz -p "$THREADS" "${RAW_DIR}/${ACC}_1.fastq" "${RAW_DIR}/${ACC}_2.fastq" 2>/dev/null || gzip -f "${RAW_DIR}/${ACC}_1.fastq" "${RAW_DIR}/${ACC}_2.fastq"
        elif [[ -f "${RAW_DIR}/${ACC}.fastq" ]]; then
            pigz -p "$THREADS" "${RAW_DIR}/${ACC}.fastq" 2>/dev/null || gzip -f "${RAW_DIR}/${ACC}.fastq"
            mv "${RAW_DIR}/${ACC}.fastq.gz" "$R1" 2>/dev/null || true
        fi
    fi

    # 2. Trim (Handles Single-End & Paired-End)
    echo "[$ACC] Trimming adapters and low-quality bases with fastp..."
    if [[ -f "$R2" && -s "$R2" ]]; then
        fastp \
            --in1 "$R1" --in2 "$R2" \
            --out1 "${TRIM_DIR}/${ACC}_1.trim.fastq.gz" \
            --out2 "${TRIM_DIR}/${ACC}_2.trim.fastq.gz" \
            --json "${QC_DIR}/${ACC}.fastp.json" \
            --html "${QC_DIR}/${ACC}.fastp.html" \
            --thread "$THREADS"
        TRIM_R1="${TRIM_DIR}/${ACC}_1.trim.fastq.gz"
        TRIM_R2="${TRIM_DIR}/${ACC}_2.trim.fastq.gz"
    else
        fastp \
            --in1 "$R1" \
            --out1 "${TRIM_DIR}/${ACC}_1.trim.fastq.gz" \
            --json "${QC_DIR}/${ACC}.fastp.json" \
            --html "${QC_DIR}/${ACC}.fastp.html" \
            --thread "$THREADS"
        TRIM_R1="${TRIM_DIR}/${ACC}_1.trim.fastq.gz"
        TRIM_R2=""
    fi

    # 3. Align (Handles Single-End & Paired-End)
    echo "[$ACC] Aligning reads with BWA-MEM..."
    BAM="${ALIGN_DIR}/${ACC}.sorted.bam"

    if [[ -n "$TRIM_R2" && -f "$TRIM_R2" ]]; then
        bwa mem -t "$THREADS" "$REF" "$TRIM_R1" "$TRIM_R2" | \
            samtools view -b | \
            samtools sort -@ "$THREADS" -o "$BAM" -
    else
        bwa mem -t "$THREADS" "$REF" "$TRIM_R1" | \
            samtools view -b | \
            samtools sort -@ "$THREADS" -o "$BAM" -
    fi
    samtools index "$BAM"

    samtools depth -a "$BAM" > "${QC_DIR}/${ACC}.depth.txt"
    samtools flagstat "$BAM" > "${QC_DIR}/${ACC}.flagstat.txt" 

    # 4. Primer trimming & Consensus
    SCHEME=$(get_primer_scheme "$ACC")
    PRIMER_BED="primer_schemes/${SCHEME}.bed" 
    WORKING_BAM="$BAM"

    if [[ "$SCHEME" != "shotgun" && "$SCHEME" != "none" && -n "$SCHEME" ]]; then
        if [[ ! -f "$PRIMER_BED" ]]; then
            echo "ERROR: Primer scheme file $PRIMER_BED not found for sample $ACC!" >&2
            exit 1
        fi

        echo "[$ACC] Trimming primers using scheme: $SCHEME..."
        TRIMMED_BAM="${ALIGN_DIR}/${ACC}.ptrim.sorted.bam"
        TRIMMED_PREFIX="${ALIGN_DIR}/${ACC}.ptrim"

        ivar trim -b "$PRIMER_BED" -i "$BAM" -p "$TRIMMED_PREFIX" -e
        samtools sort -@ "$THREADS" "${TRIMMED_PREFIX}.bam" -o "$TRIMMED_BAM"
        samtools index "$TRIMMED_BAM"

        rm -f "${ALIGN_DIR}/${ACC}.ptrim.bam"
        WORKING_BAM="$TRIMMED_BAM"
    else
        echo "[$ACC] Skipping primer trimming (Scheme: ${SCHEME:-none})..."
    fi

    echo "[$ACC] Calling consensus sequence with iVar..."
    samtools mpileup -A -d 100000 -Q 20 -a -f "$REF" "$WORKING_BAM" | \
        ivar consensus -p "${CONSENSUS_DIR}/${ACC}.consensus" -q 20 -t 0.6 -m 10 -n N

    mv "${CONSENSUS_DIR}/${ACC}.consensus.fa" "${CONSENSUS_DIR}/${ACC}.consensus.fasta" 2>/dev/null || true
    CONSENSUS_FA="${CONSENSUS_DIR}/${ACC}.consensus.fasta"

    # 5. Lineage
    echo "[$ACC] Running Nextclade..."
    nextclade run \
        --input-dataset "ref/nextclade_sars-cov-2" \
        --output-tsv "${LINEAGE_DIR}/${ACC}.nextclade.tsv" \
        "$CONSENSUS_FA"

    echo "[$ACC] Running Pangolin..."
    pangolin "$CONSENSUS_FA" --analysis-mode usher --outdir "$LINEAGE_DIR"
    
    if [[ -f "${LINEAGE_DIR}/lineage_report.csv" ]]; then
        mv "${LINEAGE_DIR}/lineage_report.csv" "${LINEAGE_DIR}/${ACC}.pangolin.csv"
    fi

    echo "[$ACC] Lineage assignment completed."

done < "$SAMPLES_FILE"

echo "Rolling up QC Outputs with MultiQC..."
multiqc "$QC_DIR" "$ALIGN_DIR" -o "$QC_DIR"

echo "Pipeline complete."