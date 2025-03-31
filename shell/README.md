# SvABA Somatic SV Calling Pipeline

## 🧬 Overview

This script performs somatic structural variant (SV) calling using **SvABA** in a Singularity environment.  
It runs **tumor-normal paired analysis** for a patient, limited to **amplicon breakpoint regions** defined in a BED file.

---

## 📦 Prerequisites

### 1. Directory Structure

Extract the provided dataset or prepare the following directory structure:
```
./glass_ms_bed
./reference
./singularity_def
svaba_pair_singularity.sh
```

### 2. Build the Singularity image

Get into the ./singularity_def directory:
& Run the following command to build the Singularity image
```
cd ./singularity_def
singularity build --fakeroot svaba_v1.2.0.sif svaba_v1.2.0.def
cd ..
```

### 3. Run the BAM slicing script
Execute the script with **the absolute path of the TM and NM BAM file**:
```
bash ./svaba_pair_singularity.sh /absolute/path/to/TM_BAM_file /absolute/path/to/NM_BAM_file
``` 

### 4. Check the output
The output files will be stored in the following directory:
./results/SvABA/${PATIENT_ID}/${TM_OBJECT_ID}__${NM_OBJECT_ID}
