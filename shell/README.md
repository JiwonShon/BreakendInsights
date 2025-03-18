# Slice BAM 

## Overview

This script slices a given BAM file into smaller regions defined in a BED file using a Singularity environment.

## Prerequisites

### 1. Extract the provided ZIP file

Extracting the ZIP file will create the following directories and files:
```
./bed
./reference
./singularity_def
slice_bam_pair_singularity.sh
```

### 2. Build the Singularity image

Get into the `./singularity_def` directory:
& Run the following command to build the Singularity image

```
cd ./singularity_def
singularity build --fakeroot bamtools_v1.3.sif bamtools_v1.3.def
cd ..
```

### 3. Run the BAM slicing script
Execute the script with the absolute path of the BAM file:
```
bash ./slice_bam_pair_singularity.sh /absolute/path/to/BAM_file
```

### 4. Check the output
The output files will be stored in the following directory:
```
./results/sliced_BAM/${PATIENT_ID}/${OBJECT_ID}
```
