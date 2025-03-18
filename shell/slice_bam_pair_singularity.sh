#!/bin/bash

START_TIME=$SECONDS

# =================
# Input arguments
# =================
BAM_FILE=$(realpath "$1")

PATIENT_ID=$(basename "$BAM_FILE" | cut -d'-' -f1,2,3)                              
BED_FILE="./bed/${PATIENT_ID}_ms_BNDregion.bed"                                     
OBJECT_ID=$(basename "$BAM_FILE" | cut -d'.' -f1)

# =================
# Change the following paths according to your setup
# =================
SINGULARITY_IMAGE="./singularity_def/bamtools_v1.3.sif"             
REFERENCE_FILE="./reference/human_g1k_v37_decoy.fasta"                       
BIND_PATHS="$(pwd):$(pwd)"                                           

# =================
# Check files
# =================
if [ ! -f "$BAM_FILE" ]; then
    echo "❌ ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi
echo "✅ BAM file loaded: $BAM_FILE"

if [ ! -f "$BED_FILE" ]; then
    echo "❌ ERROR: BED file not found: $BED_FILE"
    exit 1
fi
echo "✅ BED file loaded: $BED_FILE"

if [ ! -f "$REFERENCE_FILE" ]; then
    echo "❌ ERROR: Reference genome file not found: $REFERENCE_FILE"
    exit 1
fi
echo "✅ Reference genome loaded: $REFERENCE_FILE"


WORKDIR="$(pwd)/results/sliced_BAM/${PATIENT_ID}/${OBJECT_ID}"
mkdir -p "$WORKDIR"
echo "📁 Working directory: $WORKDIR"


# =================
# Work process
# =================
# --- Step 1: Generate BAM Header ---
echo "🔹 Generating BAM header..."
singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
    samtools view -H "$BAM_FILE" > "$WORKDIR/$(basename "$BAM_FILE").header.txt"

# --- Step 2: Generate BAM Index ---
echo "🔹 Generating BAM index..."
singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
    samtools index "$BAM_FILE"

# --- Step 3: Run BAM Slicing ---
echo "🚀 Running BAM slicing..."
singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
    java -Xmx50G -cp /opt/hmftools/bamtools/bam-tools_v1.3.jar \
    com.hartwig.hmftools.bamtools.slice.RegionSlicer \
    -bam_file "$BAM_FILE" \
    -output_dir "$WORKDIR" \
    -output_prefix "${OBJECT_ID}_BNDregion" \
    -ref_genome "$REFERENCE_FILE" \
    -regions_file "$BED_FILE" \
    -write_bam \
    -log_level INFO

# --- Step 4: Finalize ---
echo "0" > "$WORKDIR/finish_flag.txt"
ls -al -R "$WORKDIR" > "$WORKDIR/env.txt"

# --- Step 5: Compute execution time ---
END_TIME=$SECONDS
ELAPSED_TIME=$((END_TIME - START_TIME))
ELAPSED_HOUR=$((ELAPSED_TIME / 3600))
ELAPSED_MIN=$(((ELAPSED_TIME % 3600) / 60))
ELAPSED_SEC=$((ELAPSED_TIME % 60))

echo "✅ Processing completed successfully!"
echo "⏳ Total execution time: ${ELAPSED_HOUR} hr ${ELAPSED_MIN} min ${ELAPSED_SEC} sec."
