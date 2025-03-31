#!/bin/bash

START_TIME=$SECONDS

# =================
# Input arguments
# =================
TM_BAM_FILE=$(realpath "$1")
NM_BAM_FILE=$(realpath "$2")

PATIENT_ID=$(basename "$TM_BAM_FILE" | cut -d'-' -f1,2,3)                              
TM_OBJECT_ID=$(basename "$TM_BAM_FILE" | cut -d'.' -f1)
NM_OBJECT_ID=$(basename "$NM_BAM_FILE" | cut -d'.' -f1)

BAM_DIR=$(dirname "$TM_BAM_FILE")
BED_FILE="./glass_ms_bed/${PATIENT_ID}_ms_BNDregion.bed"
# =================
# Change the following paths according to your setup
# =================
SINGULARITY_IMAGE="./singularity_def/svaba_v1.2.0.sif"             
REFERENCE_FILE="./reference/human_g1k_v37_decoy.fasta"   
svaba_dbsnp_vcf="./reference/dbsnp_138.b37.vcf"
BIND_PATHS="$(pwd):$(pwd),${BAM_DIR}:${BAM_DIR}"                                   

# =================
# Check files
# =================
if [ ! -f "$BED_FILE" ]; then
    echo "❌ ERROR:BED file not found: $BED_FILE"
    exit 1
fi
echo "✅ BED file loaded: $BED_FILE"
if [ ! -f "$TM_BAM_FILE" ]; then
    echo "❌ ERROR:TM BAM file not found: $TM_BAM_FILE"
    exit 1
fi
echo "✅ TM BAM file loaded: $TM_BAM_FILE"
if [ ! -f "$NM_BAM_FILE" ]; then
    echo "❌ ERROR: NM BAM file not found: $NM_BAM_FILE"
    exit 1
fi
echo "✅ NM BAM file loaded: $NM_BAM_FILE"
if [ ! -f "$REFERENCE_FILE" ]; then
    echo "❌ ERROR: Reference genome file not found: $REFERENCE_FILE"
    exit 1
fi
echo "✅ Reference genome loaded: $REFERENCE_FILE"
if [ ! -f "$svaba_dbsnp_vcf" ]; then
    echo "❌ ERROR: Reference genome file not found: $REFERENCE_FILE"
    exit 1
fi
echo "✅ svaba_dbsnp_vcf loaded: $svaba_dbsnp_vcf"


WORKDIR="$(pwd)/results/SvABA/${PATIENT_ID}/${TM_OBJECT_ID}__${NM_OBJECT_ID}"
mkdir -p "$WORKDIR/bam"
echo "📁 Working directory: $WORKDIR"


# =================
# Work process
# =================
# --- Step 1: Generate BAM Header ---
echo "🔹 Generating BAM header..."
singularity exec "$SINGULARITY_IMAGE" which samtools

singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
    samtools view -H "$TM_BAM_FILE" > "$WORKDIR/bam/$(basename "$TM_BAM_FILE").header.txt" 
singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
    samtools view -H "$NM_BAM_FILE" > "$WORKDIR/bam/$(basename "$NM_BAM_FILE").header.txt" 

# --- Step 2: Generate BAM Index ---
# echo "🔹 Generating BAM index..."
# singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
#     samtools index "$BAM_FILE"

# --- Step 3: Run SvABA  ---
echo "🚀 Running SvABA including Amplicon interval region..."
singularity exec --bind "$BIND_PATHS" "$SINGULARITY_IMAGE" \
    svaba run \
    -t "$TM_BAM_FILE" \
    -n "$NM_BAM_FILE" \
    -D "$svaba_dbsnp_vcf" \
    -a "$WORKDIR/${TM_OBJECT_ID}__${NM_OBJECT_ID}" \
    -G "$REFERENCE_FILE" \
    -k "$BED_FILE" && 
    ls -al -R ./ >> env.txt

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
