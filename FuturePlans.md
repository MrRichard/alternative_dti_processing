# Future Plans: Batch Processing Pipeline

This document outlines a plan for extending the single-subject DTI pipeline to support
batch processing of large cohorts, with automated job scheduling and quality control.

## Overview

The current `process_dti.sh` script is designed for single-subject processing via SLURM `srun`.
A batch processing pipeline would enable:

- Automated processing of hundreds of subjects without manual intervention
- Parallel processing of multiple subjects to reduce total turnaround time
- Automated quality control and failure notification
- Centralized logging and monitoring

---

## Pipeline Architecture

### High-Level Flow

```
1. DICOM → NIfTI Conversion (dcm2niix)
   └─── Input: DICOM directory (per subject)
   └─── Output: NIfTI files + JSON sidecars

2. File Identification (SeriesDescription parsing)
   └─── Scan JSON files for SeriesDescription fields
   └─── Categorize files: AP_DWI, PA_DWI, T1w, etc.

3. Skull-Stripping (deepbet.sif for T1w)
   └─── Run deep learning brain extraction on T1w image
   └─── Output: brain-masked T1w (for potential registration)

4. SLURM Job Submission
   └─── Build `srun` command with correct inputs
   └─── Submit to appropriate partition (defq, longq, etc.)
   └─── Capture job ID for monitoring

5. Parallel Batch Processing
   └─── Run N subjects simultaneously (configurable)
   └─── Wait for all jobs to complete or fail

6. Quality Control
   └─── Automated QC metrics (SNR, FA, brain coverage)
   └─── Visual QC (generate QC snapshots for review)
   └─── Report generation (success/failure per subject)

7. Group-Level Data Collection
   └─── Aggregate normalized DTI maps
   └─── Prepare for group statistics (TBSS, ROI, voxel-wise)
```

---

## Implementation Details

### 1. DICOM to NIfTI Conversion (dcm2niix)

**Tool**: `dcm2niix` (GitHub: https://github.com/richardtan/dcm2niix)

**Input**: DICOM directory (single session)
**Output**: 
- `sub-001_AP.nii.gz` (AP DWI)
- `sub-001_AP.json`
- `sub-001_AP.bvec`, `sub-001_AP.bval`
- `sub-001_PA.nii.gz` (PA DWI)
- `sub-001_T1w.nii.gz` (T1-weighted anatomical)

**Command**:
```bash
dcm2niix \
  -b y \          # Export bvec/bval
  -z y \          # Compress to .nii.gz
  -f %s_%d \      # Filename format
  -o /output_dir  \
  /input_dicom_dir
```

**Key considerations**:
- Ensure consistent naming conventions (e.g., `sub-XXXX_YY_...`)
- Handle multi-session studies by organizing per-session directories
- Preserve DICOM headers if re-conversion is needed

### 2. File Identification via SeriesDescription

**Approach**: Parse JSON sidecar files for `SeriesDescription` field.

**Example JSON**:
```json
{
  "SeriesDescription": "NHP_DTI_30dir_iPAT4_AP",
  "ProtocolName": "DTI_30dir_iPAT4_AP",
  "AcquisitionType": "SERIES",
  "TotalReadoutTime": 0.06
}
```

**Processing logic**:
```bash
# Extract SeriesDescription from JSON
SERIES_DESC=$(python3 -c "import json; print(json.load(open('file.json'))['SeriesDescription'])")

# Categorize based on keywords
if [[ "$SERIES_DESC" == *"_AP_"* ]]; then
    # AP DWI
elif [[ "$SERIES_DESC" == *"_PA_"* ]]; then
    # PA DWI (for topup)
elif [[ "$SERIES_DESC" == *"T1w"* ]] || [[ "$SERIES_DESC" == *"MPRAGE"* ]]; then
    # T1-weighted anatomical
fi
```

**Alternative**: Use BIDS conventions if available (`_acq-ap`, `_acq-pa`, `_T1w`).

### 3. Skull-Stripping with deepbet.sif

**Tool**: [NHP-BrainExtraction](https://github.com/HumanBrainED/NHP-BrainExtraction)

**Description**: Deep learning-based brain extraction for non-human primate MRI.
Pre-trained models for macaque T1w images.

**Command**:
```bash
singularity exec deepbet.sif \
    deepbet \
    -i /input/T1w.nii.gz \
    -o /output/T1w_brain.nii.gz \
    -m /output/T1w_brain_mask.nii.gz
```

**Rationale**: 
- More robust than FSL BET for NHP T1w images
- Handles variable resolution and contrast
- Can be used for T1w→DWI registration (optional: use T1w brain mask to improve DWI mask)

**Integration**:
- Run before DTI pipeline (optional step)
- Use T1w brain mask as `ECC_MASK` for `process_dti.sh -c` (optional)
- Can improve brain extraction quality on difficult NHP datasets

### 4. SLURM Job Submission

**Approach**: Generate SLURM batch scripts dynamically.

**Example job script** (`job_<subject>.sh`):
```bash
#!/bin/bash
#SBATCH --job-name=dti_<subject>
#SBATCH --partition=defq
#SBATCH --account=ansir-users
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --output=/logs/dti_<subject>_%j.log
#SBATCH --error=/logs/dti_<subject>_%j.err

module load apptainer
srun apptainer exec \
  -B /data/dicom:/input \
  -B /data/output:/output \
  -B /templates:/templates \
  /path/to/alternative_dti_processing.sif \
  process_dti.sh \
  -a /input/${SUBJECT}_AP.nii.gz \
  -v /input/${SUBJECT}_AP.bvec \
  -b /input/${SUBJECT}_AP.bval \
  -p /input/${SUBJECT}_PA.nii.gz \
  -o /output/${SUBJECT}/ \
  -S nhp \
  -f /templates/DTIPaxFA.nii.gz
```

**SLURM array jobs** (for parallel processing):
```bash
# Submit 10 subjects in parallel
#SBATCH --array=1-10
srun apptainer exec ... process_dti.sh ... -o /output/${SLURM_ARRAY_TASK_ID}/
```

### 5. Parallel Batch Processing

**Parallelization strategy**:
- Submit multiple jobs simultaneously (configurable: 5, 10, 20 subjects)
- Use SLURM job dependencies to wait for completion
- Implement retry logic for failed jobs (2-3 retries with exponential backoff)

**Example bash loop**:
```bash
for subject in $(cat subject_list.txt); do
    ./generate_job.sh "$subject" | sbatch
    # Limit concurrent jobs to 10
    while [[ $(squeue -u $USER -h | wc -l) -ge 10 ]]; do
        sleep 30
    done
done
```

**Monitoring**:
- Track job status: `squeue -u $USER`
- Check logs: `/logs/dti_*.log`
- Send notifications on failure: `sbatch --mail-type=END,FAIL ...`

### 6. Quality Control

**Automated QC metrics**:

| Metric | Description | Threshold | Tool |
|--------|-------------|-----------|------|
| SNR | Signal-to-noise ratio in b=0 | >10 | `mrstats` |
| FA_mean | Mean fractional anisotropy | >0.2 | `mrstats` |
| brain_coverage | % brain voxels in mask | >90% | `fslmaths` |
| FA_std | FA standard deviation | <0.15 | `mrstats` |
| bias_variance | Bias field variation | <0.1 | `dwibiascorrect` output |

**QC script**:
```bash
#!/bin/bash
# qc_dti.sh - Run automated QC on DTI output

SUBJECT=$1
OUTPUT_DIR="/output/$SUBJECT"

# Check if outputs exist
for file in ${OUTPUT_DIR}/${SUBJECT}_FA.nii.gz \
            ${OUTPUT_DIR}/${SUBJECT}_mask.nii.gz; do
    if [[ ! -f "$file" ]]; then
        echo "FAIL: Missing output $file"
        exit 1
    fi
done

# Compute FA mean/std
FA_STATS=$(mrstats --mask=${OUTPUT_DIR}/${SUBJECT}_mask.nii.gz \
                  --output=mean, std \
                  ${OUTPUT_DIR}/${SUBJECT}_FA.nii.gz | tail -1)
FA_MEAN=$(echo $FA_STATS | awk '{print $1}')
FA_STD=$(echo $FA_STATS | awk '{print $2}')

echo "${SUBJECT} FA_mean=${FA_MEAN:.3f} FA_std=${FA_STD:.3f}"

# Check thresholds
if (( $(echo "$FA_MEAN < 0.2" | bc -l) )); then
    echo "WARN: Low mean FA"
fi
if (( $(echo "$FA_STD > 0.15" | bc -l) )); then
    echo "WARN: High FA variance"
fi
```

**Visual QC**:
- Generate PNG snapshots of FA maps (slice-by-slice montage)
- Overlay brain masks on anatomical (if available)
- Export to HTML report for manual review



---

## References

- **dcm2niix**: https://github.com/richardtan/dcm2niix
- **NHP-BrainExtraction**: https://github.com/HumanBrainED/NHP-BrainExtraction
- **MRtrix3**: https://www.mrtrix.org
- **FSL**: https://fsl.fmrib.ox.ac.uk
- **ANTs**: https://github.com/ANTsX/ANTs
- **SLURM**: https://slurm.schedmd.com

---

**Last updated**: June 26, 2026
