# NHP DTI Processing Pipeline — Student Guide

This guide walks through every step of the NHP diffusion tensor imaging (DTI) pipeline
(`process_dti_nhp.sh`), explaining what each command does, which software it comes from,
and where to find the official documentation.

---

## Table of Contents

1. [What is a Singularity Container?](#1-what-is-a-singularity-container)
2. [What is a SIF File?](#2-what-is-a-sif-file)
3. [Running the Container](#3-running-the-container)
4. [Pipeline Overview](#4-pipeline-overview)
5. [Step 1 — Convert to MIF Format](#step-1--convert-to-mif-format)
6. [Step 2 — Denoising (MP-PCA)](#step-2--denoising-mp-pca)
7. [Step 3 — Gibbs Ringing Correction](#step-3--gibbs-ringing-correction)
8. [Step 4 — Build the SE-EPI b0 Pair](#step-4--build-the-se-epi-b0-pair)
9. [Step 5 — Distortion and Eddy Current Correction](#step-5--distortion-and-eddy-current-correction)
10. [Step 6 — Brain Mask via FLIRT Coregistration (NHP-specific)](#step-6--brain-mask-via-flirt-coregistration-nhp-specific)
11. [Step 7 — Bias Field Correction](#step-7--bias-field-correction)
12. [Step 8 — Extract DTI Shell](#step-8--extract-dti-shell)
13. [Step 9 — Fit the Diffusion Tensor](#step-9--fit-the-diffusion-tensor)
14. [Step 10 — Extract DTI Metrics (Native Space)](#step-10--extract-dti-metrics-native-space)
15. [Step 11 — Normalize to Template Space (ANTs SyN)](#step-11--normalize-to-template-space-ants-syn)
16. [Output Files Summary](#output-files-summary)
17. [Software Versions in This Container](#software-versions-in-this-container)

---

## 1. What is a Singularity Container?

A **container** is a self-contained software environment that packages an application
together with everything it needs to run — its operating system libraries, tools, and
configuration — into a single portable file. Think of it like a shipping container: the
same contents arrive and work identically regardless of which port (computer) they land on.

**Singularity** (now also distributed as **Apptainer**) is a container system designed
specifically for scientific computing and HPC (High Performance Computing) clusters. Unlike
Docker — which is common in software engineering — Singularity was built to:

- Run without root/administrator privileges on shared cluster systems
- Integrate seamlessly with the host filesystem via **bind mounts** (making your data
  accessible inside the container without copying it)
- Maintain reproducibility: the exact same versions of FSL, MRtrix3, and ANTs run on
  every machine

When a lab publishes a Singularity container for a pipeline, you can be confident that
the analysis environment is identical to the one used to develop and validate it.

---

## 2. What is a SIF File?

A **SIF file** (Singularity Image Format, `.sif`) is the compiled, single-file form of a
Singularity container. It is built from a **definition file** (`.def`) — a plain-text
recipe that describes what to install and configure.

The workflow is:

```
Singularity_nhp.def  →  (build)  →  Singularity_nhp.sif
      (recipe)                          (executable image)
```

Once built, the `.sif` file is portable. Copy it to any Linux machine with Singularity
installed and it will run identically. You do not need to install FSL, MRtrix3, or ANTs
separately.

**Build command** (requires root or `--fakeroot` on your cluster):
```bash
singularity build Singularity_nhp.sif Singularity_nhp.def
```

---

## 3. Running the Container

The NHP container **does not run the pipeline automatically**. You must explicitly call the
script so you can control your bind mounts and arguments.

### Bind mounts (`-B`)

Your data files live on the host machine. Inside the container, the filesystem is isolated.
The `-B /host/path:/container/path` flag creates a "bridge" so the container can see your
files:

```
-B /data/nhp_study/sub01:/data
```

This means `/data/nhp_study/sub01/scan.nii.gz` on your host appears as `/data/scan.nii.gz`
inside the container.

### Example invocation

```bash
singularity run \
    -B /data/nhp_study/sub01:/data \
    -B /data/nhp_study/sub01/output:/out \
    Singularity_nhp.sif \
    /app/process_dti_nhp.sh \
        -a /data/AP_DWI.nii.gz \
        -v /data/AP.bvec \
        -b /data/AP.bval \
        -p /data/PA_DWI.nii.gz \
        -t /data/T1w.nii.gz \
        -m /data/T1w_brain_mask.nii.gz \
        -f /data/NMT_FA_template.nii.gz \
        -o /out/sub01
```

To get an interactive shell inside the container (useful for debugging):
```bash
singularity shell -B /data/nhp_study:/data Singularity_nhp.sif
```

---

## 4. Pipeline Overview

```
Raw AP DWI + PA b0  ──►  [1] mrconvert     Convert to .mif
                     ──►  [2] dwidenoise    Thermal noise removal
                     ──►  [3] mrdegibbs     Gibbs ringing removal
                     ──►  [4] dwiextract    Build b0 pair for topup
                          mrmath / mrcat
                     ──►  [5] dwifslpreproc  Distortion + eddy correction
                     ──►  [6] FLIRT          T1w mask → DWI space  ← NHP-specific
                     ──►  [7] dwibiascorrect  Bias field correction
                     ──►  [8] dwiextract     Isolate DTI shell (b=0 + b=1000)
                     ──►  [9] dwi2tensor     Fit diffusion tensor
                     ──► [10] tensor2metric  FA / MD / AD / RD (native space)
                     ──► [11] ANTs SyN       Warp to FA template  ← NHP-specific
```

The pipeline processes **multi-shell DWI** data (b = 0, 500, 1000, 2000 s/mm²) but fits
the diffusion tensor using only the **b=0 and b=1000** shells — the standard clinical DTI
approach. The full pre-processed multi-shell dataset is also saved for future advanced
modelling (NODDI, HARDI, etc.).

---

## Step 1 — Convert to MIF Format

**Software:** MRtrix3
**Documentation:** [mrconvert](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrconvert.html)

```bash
mrconvert AP_DWI.nii.gz -fslgrad AP.bvec AP.bval ap_raw.mif
mrconvert PA_DWI.nii.gz pa_raw.mif
```

`mrconvert` converts images between different file formats. Here it reads the AP DWI
series from NIfTI (`.nii.gz`) — the standard neuroimaging format — and converts it into
**MIF format** (`.mif`), MRtrix3's native format.

The `-fslgrad` flag attaches the FSL-format gradient table (`.bvec` = gradient directions,
`.bval` = b-values) directly into the MIF header. Keeping the gradient information embedded
in the file prevents the most common source of DTI errors: mismatched gradient files.

The PA series (used only for distortion correction) is converted separately without
gradient information because it is a b=0 series (no diffusion encoding).

**`mrinfo`** is also used throughout the pipeline to query image properties (dimensions,
shells, number of volumes) and write them to the log.
**Documentation:** [mrinfo](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrinfo.html)

---

## Step 2 — Denoising (MP-PCA)

**Software:** MRtrix3
**Documentation:** [dwidenoise](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwidenoise.html)

```bash
dwidenoise ap_raw.mif ap_denoised.mif -noise noise_map.mif
```

DWI images are inherently noisy because diffusion encoding gradients reduce the signal
substantially. `dwidenoise` uses **Marchenko-Pastur PCA (MP-PCA)** to separate true signal
from thermal noise.

The algorithm exploits the fact that a multi-volume DWI dataset is highly redundant:
across many gradient directions, the signal should lie in a low-dimensional subspace.
Random thermal noise, by contrast, contributes to all principal components equally. The
Marchenko-Pastur distribution gives the expected distribution of eigenvalues for a purely
random matrix, so any eigenvalues above that threshold are considered signal; those below
are noise and can be removed.

**Why this must be first:** Denoising relies on the spatial correlations in the raw
scanner data. Any interpolation (from motion correction, resampling, etc.) destroys those
correlations and makes the noise estimate incorrect. This step must always come before any
spatial manipulation.

The `-noise` flag saves a noise-level map, which can be inspected to verify uniform noise
across the brain.

---

## Step 3 — Gibbs Ringing Correction

**Software:** MRtrix3
**Documentation:** [mrdegibbs](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrdegibbs.html)

```bash
mrdegibbs ap_denoised.mif ap_degibbs.mif
```

**Gibbs ringing** (also called truncation artefact) appears as faint parallel "ringing"
lines near sharp intensity boundaries — most visibly at the brain surface and at
grey/white matter interfaces. It is caused by the finite bandwidth of the k-space
acquisition: the Fourier transform of a truncated signal produces these oscillations
(a consequence of the Gibbs phenomenon in mathematics).

`mrdegibbs` corrects for this using a local subvoxel-shift method. It must be applied
*after* denoising but *before* any geometric correction that would resample the voxels,
because the correction depends on the original k-space sampling pattern being intact.

---

## Step 4 — Build the SE-EPI b0 Pair

**Software:** MRtrix3
**Documentation:** [dwiextract](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwiextract.html) | [mrmath](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrmath.html) | [mrcat](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrcat.html)

```bash
dwiextract -bzero ap_degibbs.mif - | mrmath - mean -axis 3 b0_AP.mif
mrcat b0_AP.mif b0_PA.mif -axis 3 se_epi_pair.mif
```

To correct for EPI distortions in Step 5, the FSL `topup` tool needs a pair of b=0
(undiffusion-weighted) images acquired with **opposite phase-encoding directions** — one
AP and one PA. The distortions appear as mirror-image warps between the two, and topup
uses this to estimate and correct the field.

Three MRtrix3 commands work together here:

- **`dwiextract -bzero`** — extracts only the b=0 volumes from the AP series (volumes
  with no diffusion weighting). The `-` sends output to a pipe rather than a file.
- **`mrmath mean -axis 3`** — averages the extracted b=0 volumes along the fourth axis
  (time/volume axis) to produce a single high-SNR reference b=0. Averaging improves SNR
  proportional to the square root of the number of volumes combined.
- **`mrcat -axis 3`** — concatenates the AP b=0 and PA b=0 images into a single 4D
  image along the volume axis. This is the "spin-echo EPI pair" that `topup` requires.

---

## Step 5 — Distortion and Eddy Current Correction

**Software:** MRtrix3 wrapper / FSL (topup + eddy)
**Documentation:** [dwifslpreproc](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwifslpreproc.html)

```bash
dwifslpreproc ap_degibbs.mif ap_preproc.mif \
    -rpe_pair -se_epi se_epi_pair.mif \
    -pe_dir AP -readout_time 0.06 \
    -eddy_options "--slm=linear --data_is_shelled"
```

This is the most computationally expensive step (typically 30–90 minutes) and corrects
two types of geometric distortion:

1. **Susceptibility-induced EPI distortion** — EPI sequences are fast but highly
   susceptible to magnetic field inhomogeneities, which stretch and compress the image
   (especially in the phase-encode direction). FSL's `topup` estimates the field map
   from the AP/PA b=0 pair and computes the unwarped geometry.

2. **Eddy currents and subject motion** — Switching the strong diffusion gradients on and
   off rapidly induces currents in the scanner hardware, distorting the image. Head motion
   also occurs between volumes. FSL's `eddy` corrects both simultaneously using a model
   of the distortion field across time.

`dwifslpreproc` is a MRtrix3 **wrapper script** — it handles the orchestration between
MRtrix3 and FSL's topup/eddy, managing file conversions, calling the FSL tools in the
right order, and re-embedding the eddy-corrected (rotated) gradient directions back into
the MIF file.

Key flags:
- `-rpe_pair` — tells the script that the SE-EPI pair has reversed phase encoding
- `-pe_dir AP` — specifies the phase-encode direction of the AP series (verify against
  your scanner protocol)
- `--slm=linear` — eddy option for multi-shell data: models slice-to-volume movement with
  a linear model
- `--data_is_shelled` — tells eddy to expect a multi-shell (not single-shell) dataset

After this step, the fully corrected multi-shell DWI is exported to NIfTI
(`AP_BASE_ECC.nii.gz`) along with the eddy-rotated gradient files (`.bvec`, `.bval`).

---

## Step 6 — Brain Mask via FLIRT Coregistration (NHP-specific)

**Software:** FSL
**Documentation:** [FLIRT](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/registration/flirt)

> **Why not `dwi2mask`?**
> The standard MRtrix3 `dwi2mask` command uses a brain extraction algorithm tuned for
> human brain proportions and typical coil coverage. NHP brains differ substantially —
> smaller relative to the field of view, different skull morphology, different signal
> dropout patterns. `dwi2mask` frequently over- or under-segments NHP data. Instead, we
> use a high-quality T1w brain mask that was drawn or computed in the anatomical space,
> then coregister it into DWI space.

This step has four sub-steps:

### 6a — Extract mean corrected b0

```bash
dwiextract -bzero ap_preproc.mif - | mrmath - mean -axis 3 b0_corrected.mif
mrconvert b0_corrected.mif b0_corrected.nii.gz
```

The corrected b=0 image is used as the **registration target** (reference). It is in DWI
space and has high SNR, making it a reliable registration reference.

### 6b — FLIRT: coregister T1w → b0

```bash
flirt -in T1w.nii.gz -ref b0_corrected.nii.gz \
      -omat T1w_to_DWI.mat \
      -cost mutualinfo -dof 6 \
      -out T1w_in_DWI.nii.gz
```

**FLIRT** (FMRIB's Linear Image Registration Tool) finds the rigid-body transformation
that best aligns the T1w image to the b0. Key options:

- `-cost mutualinfo` — uses **Mutual Information** as the similarity metric. Because the
  T1w and b0 have very different contrast (T1w: bright CSF, dark WM; b=0: roughly
  opposite), a metric like correlation would fail. Mutual information measures statistical
  dependence between intensity distributions and works across modalities.
- `-dof 6` — limits the transform to **6 degrees of freedom** (3 translations + 3
  rotations, i.e. rigid body). This is appropriate because the brain geometry should not
  change between the T1w and DWI acquisitions — we are only correcting for head position.
- `-omat` — saves the transformation matrix (a 4×4 affine matrix in FSL format) for reuse
  in the next step.

### 6c — Apply transform to the brain mask

```bash
flirt -in T1w_brain_mask.nii.gz -ref b0_corrected.nii.gz \
      -applyxfm -init T1w_to_DWI.mat \
      -interp nearestneighbour \
      -out brain_mask_DWI.nii.gz
```

The same transformation matrix is applied to the binary brain mask.

- `-applyxfm -init` — applies a pre-computed transform without re-estimating it, ensuring
  the mask moves identically to the T1w image.
- `-interp nearestneighbour` — **critical for binary masks**. Standard trilinear
  interpolation would blur the edges of the mask, creating voxels with fractional values
  (e.g. 0.3, 0.7) that are not valid for a binary mask. Nearest-neighbour interpolation
  assigns each output voxel the value of the nearest input voxel, preserving the 0/1
  binary nature of the mask.

### 6d — Convert mask back to MIF

```bash
mrconvert brain_mask_DWI.nii.gz brain_mask.mif
```

Downstream MRtrix3 steps expect the mask in MIF format. This line also saves a copy to
the output directory as `AP_BASE_mask.nii.gz`.

---

## Step 7 — Bias Field Correction

**Software:** MRtrix3 (wrapper over FSL FAST)
**Documentation:** [dwibiascorrect](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwibiascorrect.html)

```bash
dwibiascorrect fsl ap_preproc.mif ap_biascorr.mif -mask brain_mask.mif
```

**B1 field inhomogeneity** (bias field) is a smooth, low-frequency variation in signal
intensity across the image caused by imperfect radiofrequency coil sensitivity profiles.
It does not reflect true tissue properties and, if uncorrected, artificially inflates or
deflates apparent diffusivity values across space.

`dwibiascorrect` is a MRtrix3 wrapper that applies FSL's `fast` (or optionally ANTs N4)
to estimate and correct this field. The `-mask` option confines the correction to brain
voxels, preventing the background noise from biasing the field estimate.

---

## Step 8 — Extract DTI Shell

**Software:** MRtrix3
**Documentation:** [dwiextract](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwiextract.html)

```bash
dwiextract -shells 0,1000 ap_biascorr.mif dwi_dti_shell.mif
```

The full dataset contains four shells (b = 0, 500, 1000, 2000 s/mm²). The diffusion
tensor model (used in Steps 9–10) assumes a single-exponential signal decay and is only
valid when using data from **two shells**: one unweighted (b=0) and one diffusion-weighted.

Using b=1000 s/mm² is the standard clinical DTI choice:
- b=500 is too low to provide good contrast between tissues
- b=2000 provides better sensitivity to restricted diffusion but deviates from the tensor
  model assumptions, and is better suited to HARDI or NODDI modelling

`dwiextract -shells` selects only the volumes belonging to the specified b-values, passing
through the gradient information for the kept volumes.

---

## Step 9 — Fit the Diffusion Tensor

**Software:** MRtrix3
**Documentation:** [dwi2tensor](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2tensor.html)

```bash
dwi2tensor dwi_dti_shell.mif dt.mif -mask brain_mask.mif
```

The **diffusion tensor** is a 3×3 symmetric matrix that describes the direction and
magnitude of water diffusion in each voxel. In white matter tracts, water diffuses
preferentially along the fibre direction (anisotropic diffusion), so the tensor is
elongated along the axon axis.

`dwi2tensor` performs a **weighted least-squares fit** of the tensor model to the
measured DWI signal. The signal equation is:

```
S(b,g) = S₀ · exp(-b · gᵀDg)
```

where `S₀` is the b=0 signal, `b` is the b-value, `g` is the gradient direction unit
vector, and `D` is the 3×3 diffusion tensor being estimated. At minimum, 7 volumes are
needed (6 non-collinear directions + b=0); more directions improve the fit.

The output `dt.mif` stores 6 independent tensor elements per voxel (the tensor is
symmetric, so 6 of the 9 elements are unique). The `-mask` option skips fitting outside
the brain, saving computation.

---

## Step 10 — Extract DTI Metrics (Native Space)

**Software:** MRtrix3
**Documentation:** [tensor2metric](https://mrtrix.readthedocs.io/en/latest/reference/commands/tensor2metric.html)

```bash
tensor2metric dt.mif \
    -fa  fa.mif  \
    -adc md.mif  \
    -ad  ad.mif  \
    -rd  rd.mif  \
    -mask brain_mask.mif
```

`tensor2metric` diagonalises the tensor at each voxel (finding its three eigenvalues
λ₁ ≥ λ₂ ≥ λ₃) and computes summary metrics:

| Flag | Metric | Formula | Interpretation |
|------|--------|---------|----------------|
| `-fa` | **Fractional Anisotropy** | `FA = √(3/2) · √[(λ₁−λ̄)²+(λ₂−λ̄)²+(λ₃−λ̄)²] / √(λ₁²+λ₂²+λ₃²)` | 0 = isotropic, 1 = fully anisotropic. Reflects white matter organisation. |
| `-adc` | **Mean Diffusivity (MD)** | `MD = (λ₁+λ₂+λ₃)/3` | Overall magnitude of diffusion. Elevated in oedema; reduced in acute ischaemia. |
| `-ad` | **Axial Diffusivity (AD)** | `AD = λ₁` | Diffusion along the principal axis; sensitive to axonal integrity. |
| `-rd` | **Radial Diffusivity (RD)** | `RD = (λ₂+λ₃)/2` | Diffusion perpendicular to fibres; sensitive to myelin integrity. |

Each metric is saved as a NIfTI file in the **native DWI space** of this subject.

---

## Step 11 — Normalize to Template Space (ANTs SyN)

**Software:** ANTs (Advanced Normalization Tools)
**Documentation:** [ANTs Wiki](https://github.com/ANTsX/ANTs/wiki)

> **Why normalize?**
> Each subject's brain sits in its own native space — a coordinate system defined by
> where the animal's head happened to be in the scanner. To compare FA or MD values
> across subjects (or against an atlas), all images must be brought into a **common
> template space** where equivalent anatomical locations share the same coordinates.

### 11a — Register native FA → FA template (SyN)

```bash
antsRegistrationSyN.sh \
    -d 3 \
    -f FA_template.nii.gz \
    -m sub01_FA.nii.gz \
    -o ants_ \
    -t s
```

`antsRegistrationSyN.sh` is a convenient wrapper script provided with ANTs that runs a
**two-stage registration**: first an affine (linear) alignment, then a **SyN
(Symmetric Normalization)** deformable registration.

- **Affine stage** — corrects for global differences in position, orientation, and scale
  (12 degrees of freedom).
- **SyN stage** (`-t s`) — estimates a smooth, invertible **deformation field** that
  locally deforms the moving image to match the template, accounting for individual
  anatomical shape differences. SyN is considered the gold standard for brain
  normalization because it is diffeomorphic (the deformation is smooth and invertible,
  meaning anatomy is never folded or torn).

Outputs:
- `ants_0GenericAffine.mat` — the affine transform
- `ants_1Warp.nii.gz` — the deformation field (native → template)
- `ants_1InverseWarp.nii.gz` — the inverse deformation field (template → native)

### 11b — Apply warp to all DTI metrics

```bash
antsApplyTransforms \
    -d 3 \
    -i sub01_FA.nii.gz \
    -r FA_template.nii.gz \
    -t ants_1Warp.nii.gz ants_0GenericAffine.mat \
    -o wSub01_ECC_FA.nii.gz \
    -n Linear
```

`antsApplyTransforms` applies the transforms computed in the previous step to each DTI
metric. Using the same registration transforms for all metrics ensures every map (FA, MD,
AD, RD) is in perfect spatial alignment in template space.

The transforms are applied in **right-to-left order**: affine first, then the deformation
field — this is the standard ANTs convention.

- `-n Linear` — trilinear interpolation is appropriate for continuous-valued scalar maps
  like FA and MD (in contrast to the nearest-neighbour interpolation used for the binary
  mask in Step 6).

The output files are prefixed with `w` (for "warped") to distinguish them from the
native-space maps.

---

## Output Files Summary

All outputs appear in the directory you specified with `-o`.

### Native DWI space

| File | Contents |
|------|----------|
| `<AP_BASE>_ECC.nii.gz` | Full preprocessed multi-shell DWI (all b-values) |
| `<AP_BASE>_ECC.bvec` | Eddy-rotated gradient directions |
| `<AP_BASE>_ECC.bval` | B-values |
| `<AP_BASE>_mask.nii.gz` | Brain mask in DWI space (derived from T1w via FLIRT) |
| `<AP_BASE>_FA.nii.gz` | Fractional Anisotropy |
| `<AP_BASE>_MD.nii.gz` | Mean Diffusivity |
| `<AP_BASE>_AD.nii.gz` | Axial Diffusivity |
| `<AP_BASE>_RD.nii.gz` | Radial Diffusivity |

### Template space (ANTs SyN normalized)

| File | Contents |
|------|----------|
| `w<AP_BASE>_ECC_FA.nii.gz` | FA warped to FA template |
| `w<AP_BASE>_ECC_MD.nii.gz` | MD warped to FA template |
| `w<AP_BASE>_ECC_AD.nii.gz` | AD warped to FA template |
| `w<AP_BASE>_ECC_RD.nii.gz` | RD warped to FA template |

---

## Software Versions in This Container

| Software | Version | Purpose |
|----------|---------|---------|
| **MRtrix3** | latest (from `mrtrix3/mrtrix3:latest`) | DWI processing, format conversion, tensor fitting |
| **FSL** | bundled with MRtrix3 base image | topup, eddy, FLIRT |
| **ANTs** | 2.6.5 | SyN registration, applying transforms |
| **Python 3** | Debian 12 default | JSON sidecar parsing (readout time) |

The container is built on **Debian 12 (bookworm)**, amd64. ANTs is the precompiled
Ubuntu 22.04 build, which is fully compatible with Debian 12's glibc 2.36 (Ubuntu 22.04
uses glibc 2.35; newer glibc versions are backward-compatible with older-built binaries).
