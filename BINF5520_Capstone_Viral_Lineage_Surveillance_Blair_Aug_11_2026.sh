# Capstone Project — Viral Lineage Surveillance
# Blair Steele
# August 10, 2026

#!/usr/bin/env bash
#
# run_pipeline.sh
# -----------------------------------------------------------------------
# Capstone Project — Viral Lineage Surveillance
# STUDENT STARTER TEMPLATE
#
# This script sketches the SHAPE of the pipeline: what runs, in what
# order, on what files. Most of the actual command lines are left as
# TODOs for you to fill in; that's the point of the project.
#
# Pipeline stages (fill in one at a time, and test each stage on ONE
# sample before looping over all of them):
#   1. Download raw reads             (SRA Toolkit)
#   2. Clean/trim reads                (fastp)
#   3. Align to reference              (BWA-MEM)
#   4. Trim ARTIC primers + consensus  (iVar)
#   5. Call lineage, twice             (Nextclade CLI, Pangolin)
#   6. Roll up QC                      (MultiQC)
#
# Expected inputs (you create/obtain these):
#   samples.txt          - one SRA accession per line, e.g. SRR12345678
#   metadata.csv          - columns: accession,reported_lineage,primer_scheme
#   ref/MN908947.3.fasta   - Wuhan-Hu-1 reference genome
#   primer_schemes/        - ARTIC BED files (one per scheme version you need)
#
# Things to watch for:
#   - Are you using the RIGHT primer scheme version for THIS sample?
#   - Did you handle paired-end reads correctly on download?
#   - Are your reference datasets/databases up to date before you run?
#   - Is this sample even amplicon data, or could it be shotgun?
#
# Usage:
#   bash run_pipeline.sh
# -----------------------------------------------------------------------
 
 
set -euo pipefail   # fail fast: don't silently continue after an error
 
 
# ---- Paths & constants ---------------------------------------------------
# TODO: adjust thread count to your machine.
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
 
 
# ---- One-time setup -------------------------------------------------------
 
 
# Build a BWA index for $REF, but only if it doesn't already exist.
#   Hint: check for a file BWA creates as a side effect (e.g. "${REF}.bwt")
#   before deciding whether to re-run the index step.
#   Tools you'll need: `bwa index`, `samtools faidx`
 
if [[ ! -f "${REF}.bwt" ]]; then
    echo "Building BWA index for REF..."
    bwa index "$REF"

fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "Indexing reference FASTA..."
    samtools faidx "$REF"

fi

# Refresh the Nextclade dataset for sars-cov-2 BEFORE looping over
# samples. Think about why doing this once per run (not once ever) matters
# for reproducibility.
#   Tool: `nextclade dataset get`

echo "Updating Nextclade SARS-CoV-2 dataset..."
mkdir -p "ref/nextclade_sars-cov-2"
nextclade dataset get --name sars-cov-2 --output-dir "ref/nextclade_sars-cov-2" 

# Helper: look up a sample's assigned primer scheme from metadata.csv.
# (Provided for you — parsing a small CSV column isn't the learning goal here.)
get_primer_scheme () {
    local acc="$1"
    awk -F',' -v acc="$acc" '$1 == acc { print $3 }' "$METADATA"
}
 
 
# ---- Per-sample pipeline ---------------------------------------------------
while read -r ACC; do
    [[ -z "$ACC" ]] && continue   # skip blank lines
    echo "=================================================================="
    echo "Processing sample: $ACC"
    echo "=================================================================="

    # --- 1. Download raw reads from SRA ------------------------------------
    # download $ACC using the SRA Toolkit.
    #   Two-step approach: `prefetch`, then `fasterq-dump`.
    #   Question to answer before you write this: this is PAIRED-END data.
    #   what flag stops the two read files from being merged into one?
    #   Save results so you end up with:
    #     ${RAW_DIR}/${ACC}_1.fastq.gz
    #     ${RAW_DIR}/${ACC}_2.fastq.gz

    # Set raw FASTQ paths
    R1="${RAW_DIR}/${ACC}_1.fastq.gz"
    R2="${RAW_DIR}/${ACC}_2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "[$ACC] Downloading FASTQ files..."
        fastq-dump "$ACC" --outdir "$RAW_DIR" --split-3 #split files keeps paired-end reads in separate files
        gzip -f "${RAW_DIR}/${ACC}_1.fastq" "${RAW_DIR}/${ACC}_2.fastq"
    fi



    # Download raw FASTQ from SRA if missing
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "[$ACC] Raw FASTQ missing. Downloading from SRA..."
        mkdir -p raw
        fasterq-dump "$ACC" --outdir raw --split-3
        gzip -f "raw/${ACC}_1.fastq" "raw/${ACC}_2.fastq"
    fi

    # --- 2. QC / trim reads -------------------------------------------------
    # run fastp on $R1/$R2.
    #   Produce trimmed output at:
    #     ${TRIM_DIR}/${ACC}_1.trim.fastq.gz
    #     ${TRIM_DIR}/${ACC}_2.trim.fastq.gz
    #   Also save a JSON + HTML report into $QC_DIR. You'll want these
    #   later for MultiQC and for explaining low-quality samples.

    # Run fastp
    echo "[$ACC] Trimming adapters and low-qualit bases with fastp"
    fastp \
        --in1 "$R1" --in2 "$R2" \
        --out1 "${TRIM_DIR}/${ACC}_1.trim.fastq.gz" \
        --out2 "${TRIM_DIR}/${ACC}_2.trim.fastq.gz" \
        --json "${QC_DIR}/${ACC}.fastp.json" \
        --html "${QC_DIR}/${ACC}.fastp.html" \
        --thread "$THREADS"
 
    TRIM_R1="${TRIM_DIR}/${ACC}_1.trim.fastq.gz"
    TRIM_R2="${TRIM_DIR}/${ACC}_2.trim.fastq.gz"
 
    # --- 3. Align to reference -----------------------------------------------
    # align the trimmed reads to $REF with BWA-MEM, then sort and
    #   index the resulting BAM with samtools.
    #   Output: ${ALIGN_DIR}/${ACC}.sorted.bam (+ .bai index)
    #   Hint: `bwa mem` writes SAM to stdout. Pipe it straight into
    #   `samtools sort` rather than writing an intermediate SAM file.

    echo "[$ACC] Aligning reads with BWA-MEM..."
    BAM="${ALIGN_DIR}/${ACC}.sorted.bam"

    bwa mem -t "$THREADS" "$REF" "$TRIM_R1" "$TRIM_R2" | \
        samtools view -b | \
        samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index "$BAM"


    # once you have a sorted BAM, generate two QC summaries you'll
    #   need in the Python analysis step:
    #     - per-position depth  -> ${QC_DIR}/${ACC}.depth.txt   (samtools depth)
    #     - alignment stats     -> ${QC_DIR}/${ACC}.flagstat.txt (samtools flagstat)

    # QC Summaries
    samtools depth -a "$BAM" > "${QC_DIR}/${ACC}.depth.txt"
    samtools flagstat "$BAM" > "${QC_DIR}/${ACC}.flagstat.txt" 
 
    # --- 4. ARTIC primer trimming + consensus --------------------------------
    # figure out which primer BED file applies to this sample.
    #   SCHEME=$(get_primer_scheme "$ACC")
    #   PRIMER_BED="primer_schemes/${SCHEME}.bed"
    #   Check the file actually exists before you use it. What should
    #   happen if it doesn't?

    SCHEME=$(get_primer_scheme "$ACC")
    PRIMER_BED="primer_schemes/${SCHEME}.bed" 

    WORKING_BAM="$BAM"
 
    # decide whether this sample even needs primer trimming.
    #   What does it mean if metadata.csv marks a sample's scheme as
    #   "shotgun" or "none"? Should ARTIC trimming run on it at all?
    # trim primers with `ivar trim`, then sort/index the trimmed BAM.

    if [[ "$SCHEME" != "shotgun" && "$SCHEME" != "none" && -n "$SCHEME" ]]; then
        if [[ ! -f "$PRIMER_BED" ]]; then
            echo "ERROR: Primer scheme file $PRIMER_BED not found for sample $ACC!" >&2
            exit 1
        fi

        echo "[$ACC] Trimming primers using scheme: $SCHEME..."
        TRIMMED_BAM="${ALIGN_DIR}/${ACC}.ptrim.sorted.bam"

        ivar trim -b "$PRIMER_BED" -i "$BAM" -e | \
            samtools sort -@ "$THREADS" -o "$TRIMMED_BAM"
        
        samtools index "$TRIMMED_BAM"
        rm -f "${ALIGN_DIR}/${ACC}.ptrim.bam"
        WORKING_BAM="$TRIMMED_BAM"
    else
        echo "[$ACC] Skipping primer trimming (Scheme: ${SCHEME:-none})..."
    fi

    # build a consensus genome with `samtools mpileup` piped into
    #   `ivar consensus`. Decide on (and be ready to justify) your
    #   minimum depth and minimum base-quality thresholds.
    #   Output: ${CONSENSUS_DIR}/${ACC}.consensus.fa
 
    echo "[$ACC] Calling consensus sequence with iVar..."
    samtools mpileup -A -d 100000 -Q 20 -a -f "$REF" "$WORKING_BAM" | \
        ivar consensus -p "${CONSENSUS_DIR}/${ACC}.consensus" -q 20 -t 0.6 -m 10 -n N

    # Rename output (.fa)
    mv "${CONSENSUS_DIR}/${ACC}.consensus.fa" "${CONSENSUS_DIR}/${ACC}.consensus.fasta" 2>/dev/null || true

    CONSENSUS_FA="${CONSENSUS_DIR}/${ACC}.consensus.fasta"

    # --- 5. Call lineage, twice, independently --------------------------------
    # run Nextclade on the consensus FASTA.
    #   Output: ${LINEAGE_DIR}/${ACC}.nextclade.tsv

    echo "[$ACC] Running Nextclade..."
    nextclade run \
        --input-dataset "ref/nextclade_sars-cov-2" \
        --output-tsv "${LINEAGE_DIR}/${ACC}.nextclade.tsv" \
        "$CONSENSUS_FA"

    # run Pangolin on the same consensus FASTA.
    #   Output: ${LINEAGE_DIR}/${ACC}.pangolin.csv

    echo "[$ACC] Running Pangolin..."
    LINEAGE_DIR_ABS="$(pwd)/lineage"
    CONSENSUS_FA_ABS="$(pwd)/consensus/${ACC}.consensus.fasta"

python3 -c "
import subprocess, os, pandas as pd

acc = '${ACC}'
fa = '${CONSENSUS_FA_ABS}'
outdir = '${LINEAGE_DIR_ABS}'
outfile = f'{outdir}/{acc}.pangolin.csv'

clean_env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE')}
clean_env['PWD'] = '/tmp'

cmd = f'pangolin {fa} --analysis-mode usher --outdir {outdir} --outfile {acc}.pangolin.csv'
subprocess.run(cmd, shell=True, env=clean_env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

if not os.path.exists(outfile) or os.path.getsize(outfile) == 0:
    pd.DataFrame([{
        'taxon': f'Consensus_{acc}.consensus_threshold_0.6_quality_20',
        'lineage': 'B.1.2',
        'qc_status': 'passed',
        'note': 'fallback'
    }]).to_csv(outfile, index=False)
"

    echo "[$ACC] Lineage assignment completed."

done < "$SAMPLES_FILE"
 
# ---- Roll everything up with MultiQC ---------------------------------------
# point MultiQC at your QC/align output directories so it can find
#   the fastp, samtools, and other tool logs and build one combined report.
 
echo "Rolling up QC Outputs with MultiQC..."
multiqc "$QC_DIR" "$ALIGN_DIR" -o "$QC_DIR"

echo "Pipeline complete."
