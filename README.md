# Alternative DTI Processing Pipeline

A containerized diffusion tensor imaging (DTI) preprocessing and tensor-fitting
pipeline for multi-shell AP/PA acquisitions, with **species-aware brain
masking** for both human and non-human primate (NHP / monkey) data.

The pipeline runs end-to-end inside a single [Apptainer/Singularity](https://apptainer.org)
container built on top of the official [MRtrix3](https://www.mrtrix.org) image,
and bundles **FSL**, **AFNI**, and **ANTs**.

---

## What it does

Starting from raw AP-encoded multi-shell DWI plus a reverse phase-encoded (PA)
b0 series, the pipeline produces eddy/susceptibility-corrected DWI and the four
standard DTI scalar maps (FA, MD, AD, RD) in native space, with optional
warping to a template.

| Step | Operation | Tool |
|-----:|-----------|------|
| 1 | Convert to `.mif` | MRtrix3 `mrconvert` |
| 2 | Denoise (MP-PCA) | MRtrix3 `dwidenoise` |
| 3 | Gibbs-ringing removal | MRtrix3 `mrdegibbs` |
| 4 | Build SE-EPI b0 pair | MRtrix3 `dwiextract`/`mrcat` |
| 5 | Topup + eddy correction | FSL via `dwifslpreproc` |
| 6 | **Brain mask (species-aware)** | `dwi2mask` *(human)* / AFNI `3dSkullStrip -monkey -blur_fwhm 2 -no_touchup + erosion` *(nhp)* |
| 7 | Bias-field correction | `dwibiascorrect fsl` |
| 8 | Extract b0 + DTI shell | MRtrix3 `dwiextract` |
| 9 | Tensor fit | MRtrix3 `dwi2tensor` |
| 10 | FA/MD/AD/RD maps | MRtrix3 `tensor2metric` |
| 11 | *(optional)* warp maps to template | ANTs SyN |

### Brain masking: the key difference between species

Accurate brain masking on monkey EPI data is hard — `dwi2mask` and FSL `bet`
both proved unreliable. The pipeline therefore selects the masking strategy
from the `-S` flag:

- **`-S human`** (default): MRtrix3 `dwi2mask` on the preprocessed DWI.
- **`-S nhp`**: FSL averages the b=0 volumes from the eddy-corrected series into one
  volume, AFNI `3dSkullStrip -monkey` extracts the brain (a preset tuned for
  NHP EPI contrast), with `-blur_fwhm 2` to stabilise the surface boundary on
  noisy EPI data and `-no_touchup` to prevent reclaiming non-brain edge voxels.
  A final 1-voxel erosion removes the residual non-brain halo.

This single mask then drives bias correction, tensor fitting, and metric
extraction — there is no separate `bet`/`dwi2mask` call downstream.

---

## Build the container

```bash
apptainer build alternative_dti_processing.sif Singularity.def
# (add --fakeroot if you are not root)
```

The build pulls MRtrix3 (base image), installs FSL verification, AFNI
(`linux_openmp_64`), and ANTs 2.6.5. It downloads ~1 GB of AFNI binaries plus
the ANTs release, so allow 20–30 minutes and a working network connection.

---

## Usage

### Command-line invocation (inside the container)

The script is invoked with a single `srun` command. Here is an example invocation
for a monkey subject using a single node with the default partition:

```bash
srun -p defq -A ansir-users singularity exec \
  -B ./sample_input_data/sample_nifti_data/:/input \
  -B /isilon/datalake/cam_can_meg/original/TrainingLabResources/sample_output_data/:/output \
  -B /isilon/datalake/cam_can_meg/original/TrainingLabResources/template_resources/UWDTIRhesusWMAtlas:/templates \
  alternative_dti_processing/alternative_dti_processing.sif \
  alternative_dti_processing/process_dti.sh \
  -a /input/NHP_1484_20230713_ed2d_diff_wfu_DTI_30dir_iPAT4_AP_20230713074902_17.nii \
  -v /input/NHP_1484_20230713_ed2d_diff_wfu_DTI_30dir_iPAT4_AP_20230713074902_17.bvec \
  -b /input/NHP_1484_20230713_ed2d_diff_wfu_DTI_30dir_iPAT4_AP_20230713074902_17.bval \
  -p /input/NHP_1484_20230713_ed2d_diff_wfu_DTI_30dir_iPAT4_PA_20230713074902_11.nii \
  -o /output/NHP_1484_20230713/ \
  -S nhp \
  -f /templates/DTIPaxFA.nii.gz
```

#### Input path mapping

- **Input data**: Map your DICOM or NIfTI input directory to `/input` inside the container.
- **Output data**: Map your output directory to `/output` inside the container.
- **Template data**: Map template resources (FA template, atlases) to `/templates` inside the container.

#### Required arguments

| Flag | Description |
|------|-------------|
| `-a FILE` | AP-encoded DWI series (NIfTI) |
| `-v FILE` | AP bvec file |
| `-b FILE` | AP bval file |
| `-p FILE` | PA-encoded b0 series for topup correction (NIfTI) |
| `-o DIR`  | Output directory (created if it does not exist) |

#### Optional arguments

| Flag | Description |
|------|-------------|
| `-S SPEC` | Species / masking strategy: `human` (default) or `nhp` (AFNI `-monkey` masking) |
| `-f FILE` | FA template (NIfTI) → warp DTI maps to template via ANTs SyN |
| `-r FLOAT` | Total EPI readout time in seconds [default: 0.06, auto-read from JSON sidecar] |
| `-s INT` | B-value shell for tensor fitting [default: 1000; use 2000 for high-b DTI] |
| `-c FILE` | Pre-made brain mask (NIfTI) → resume from bias correction onward |
| `-k` | Keep intermediate files in `<outdir>/tmp/` (default: delete) |
| `-h` | Show this help message |

#### Example invocations

##### Human data (default masking)

```bash
srun -p defq -A ansir-users singularity exec \
  -B /data/study:/input \
  -B /data/out:/output \
  alternative_dti_processing.sif \
  alternative_dti_processing/process_dti.sh \
  -a /input/sub-01_AP.nii.gz \
  -v /input/sub-01_AP.bvec \
  -b /input/sub-01_AP.bval \
  -p /input/sub-01_PA.nii.gz \
  -o /output/sub-01
```

##### Monkey / NHP data (AFNI `-monkey` masking)

```bash
srun -p defq -A ansir-users singularity exec \
  -B /data/nhp_study:/input \
  -B /data/out:/output \
  -B /templates:/templates \
  alternative_dti_processing.sif \
  alternative_dti_processing/process_dti.sh \
  -a /input/NHP_001_AP.nii.gz \
  -v /input/NHP_001_AP.bvec \
  -b /input/NHP_001_AP.bval \
  -p /input/NHP_001_PA.nii.gz \
  -o /output/NHP_001 \
  -S nhp \
  -f /templates/DTIPaxFA.nii.gz
```

##### Resuming with a hand-edited mask

If an auto-generated mask needs manual correction, edit the mask file and rerun with `-c`:

```bash
srun -p defq -A ansir-users singularity exec \
  -B /data/nhp_study:/input \
  -B /data/out:/output \
  alternative_dti_processing.sif \
  alternative_dti_processing/process_dti.sh \
  -a /input/NHP_001_AP.nii.gz \
  -v /input/NHP_001_AP.bvec \
  -b /input/NHP_001_AP.bval \
  -p /input/NHP_001_PA.nii.gz \
  -o /output/NHP_001 \
  -c /output/NHP_001/NHP_001_AP_mask_edited.nii.gz \
  -S nhp
```

This re-imports the existing ECC output and your mask and resumes from bias
correction onward, skipping topup/eddy and auto-masking.

---

## Outputs (`<outdir>/`)

| File | Description |
|------|-------------|
| `<base>_ECC.{nii.gz,bvec,bval}` | Preprocessed full multi-shell DWI (eddy-rotated gradients) |
| `<base>_mask.nii.gz` | Binary brain mask used throughout |
| `<base>_{FA,MD,AD,RD}.nii.gz` | DTI maps in native DWI space |
| `w<base>_ECC_{FA,MD,AD,RD}.nii.gz` | *(with `-f`)* DTI maps in template space |
| `pipeline.log` | Full processing log |

---

## Pipeline steps (detailed overview)

1. **Input conversion**: Raw NIfTI images are converted to MRtrix3 `.mif` format.
2. **Denoising**: Marchenko–Pastur PCA denoising (`dwidenoise`) applied before any spatial operation.
3. **Gibbs ringing removal**: Local subvoxel-shifts method (`mrdegibbs`) suppresses ringing artifacts.
4. **SE-EPI b0 pair**: Mean b=0 from AP and PA series are concatenated for topup.
5. **Distortion & eddy correction**: FSL `topup` + `eddy` via `dwifslpreproc` with `-rpe_pair`.
6. **Brain masking**:
   - Human: MRtrix3 `dwi2mask` on the preprocessed DWI.
   - NHP: FSL mean b=0 + AFNI `3dSkullStrip -monkey -blur_fwhm 2 -no_touchup + 1-voxel erosion`.
7. **Bias-field correction**: FSL FAST (`dwibiascorrect fsl`) within the brain mask.
8. **Shell selection**: `dwiextract` isolates b=0 + chosen DTI shell (default b=1000).
9. **Tensor fitting**: MRtrix3 `dwi2tensor` with iteratively reweighted linear least squares.
10. **Metric extraction**: `tensor2metric` produces FA, MD, AD, RD maps.
11. **Template normalization** *(optional, `-f`)*: ANTs SyN registers native FA → template; composite transform warps all four metrics.

---

## Multi-shell data notes

The pipeline supports multi-shell acquisitions (b = 0, 500, 1000, 2000 s/mm²):

- **Tensor fitting**: Uses b=0 + b=1000 by default (`-s 2000` for high-b DTI).
- **Full series**: The preprocessed ECC series is retained for advanced models:
  - NODDI, SMT, MSMT-CSD
  - Q-ball or CSD on individual shells
  - Hybrid models combining multiple shells

---

## Phase-encode direction

Step 5 assumes AP = `j-` (`-pe_dir AP`). Verify your phase-encode direction with
`mrinfo` or your scanner protocol sheet. If your AP series is actually
posterior→anterior, change `-pe_dir PA`.

---

## Readout time

The pipeline reads `TotalReadoutTime` from the AP JSON sidecar automatically when
present; override with `-r` if your JSON is missing or the field is misconfigured.

Derivation: `EffectiveEchoSpacing_s * (PhaseEncodeLines - 1)`

Check your DICOM header or scanner protocol sheet for the correct value.

---

## Re-running after container rebuild

If you rebuild the SIF file to update the pipeline, simply re-run the same
`srun` command — no changes to the invocation are required.

To force recomputation of the slow topup/eddy steps, delete the
`<outdir>/*_ECC.*` files or point `-o` at a fresh directory.

---

## Documentation

- [`methods_nhp.md`](methods_nhp.md) — methods-section write-up of the NHP pipeline.
- [`NHP_DTI_Pipeline_Guide.md`](NHP_DTI_Pipeline_Guide.md) — step-by-step student guide.
- [`FuturePlans.md`](FuturePlans.md) — outline for batch-processing infrastructure.
