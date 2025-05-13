#!/usr/bin/env bash


# Modify paths here ###########################################################################

ncores=4
FASTQ_LIST=/path/to/samples.list # Tab separated file with 2 columns: first column is path to the FastQ file, (either R1 or R2) and second column is the sample name
PATH_BASE=/path/to/working/dir
PATH_DB=/path/to/metaphlan/db # database will be downloaded here
REFERENCES_ALL=/path/to/your/reference/genomic/fasta/




###########################################################################
###########################################################################

PATH_SAM=${PATH_BASE}/sam/
mkdir ${PATH_BASE}
mkdir ${PATH_BASE}/sam
mkdir ${PATH_BASE}/bowtie2
mkdir ${PATH_BASE}/profiles
mkdir ${PATH_BASE}/db_markers
mkdir ${PATH_BASE}/strainphlan_output_tuned
mkdir ${PATH_BASE}/strainphlan_output_tuned/consensus_markers


MIN_BASE_COVERAGE=1 # Default 1
MIN_BASE_QUALITY=30 # Default 30
BREADTH_THRESHOLD=60 # Default 80
MIN_READS_ALIGNING=4 # Default 8

MARKER_IN_N_SAMPLES=0 # Default 80
SAMPLE_WITH_N_MARKERS=20 # Default 80


# Run Metaphlan ###################################################

counter=1
end=0
start=0

while IFS=$'\t' read -r f sample_id
do
    echo "Time taken for previous iteration: $(( end - start )) seconds"
    echo "Running iteration $counter"
    counter=$(( counter + 1 ))
    start=$(date +%s)

    echo "Processing sample: $sample_id (file: $f)"

    metaphlan "$f" \
      -x mpa_vFeb24_CDIFF_CHOCOPhlAnSGB_20240910 \
      --bowtie2db "$PATH_DB" \
      -s "${PATH_BASE}/sam/${sample_id}.sam.bz2" \
      --bowtie2out "${PATH_BASE}/bowtie2/${sample_id}.bowtie2.bz2" \
      --nproc "$ncores" \
      --input_type fastq \
      --unclassified_estimation \
      -o "${PATH_BASE}/profiles/${sample_id}.txt"

    end=$(date +%s)
done < "$FASTQ_LIST"


## Extract Markers ###################################################

extract_markers.py \
    -c t__GGB4453 \
    -d ${PATH_DB}/mpa_vFeb24_CDIFF_CHOCOPhlAnSGB_20240910.pkl \
    -o ${PATH_BASE}/db_markers

## Get sample markers ###################################################

sample2markers.py \
    -i ${PATH_BASE}/sam/*.sam.bz2 \
    -o ${PATH_BASE}/strainphlan_output_tuned/consensus_markers \
    -n ${ncores} \
    -d ${PATH_DB}/mpa_vFeb24_CDIFF_CHOCOPhlAnSGB_20240910.pkl \
    --min_reads_aligning ${MIN_READS_ALIGNING} \
    --min_base_coverage ${MIN_BASE_COVERAGE} \
    --min_base_quality ${MIN_BASE_QUALITY} \
    --breadth_threshold ${BREADTH_THRESHOLD}

## Run Strainphlan ###################################################

strainphlan \
    -s ${PATH_BASE}/strainphlan_output_tuned/consensus_markers/*.pkl \
    -d ${PATH_DB}/mpa_vFeb24_CDIFF_CHOCOPhlAnSGB_20240910.pkl \
    -r ${REFERENCES_ALL} \
    -o ${PATH_BASE}/strainphlan_output_br${BREADTH_THRESHOLD}_mra${MIN_READS_ALIGNING} \
    -n ${ncores} \
    -c t__GGB4453 \
    --mutation_rates \
    --sample_with_n_markers ${SAMPLE_WITH_N_MARKERS} \
    --marker_in_n_samples ${MARKER_IN_N_SAMPLES} \
    --breadth_thres ${BREADTH_THRESHOLD}
