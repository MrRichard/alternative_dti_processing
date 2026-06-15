# Alternative DTI Processing

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
| 6 | **Brain mask (species-aware)** | `dwi2mask` *(human)* / AFNI `3dSkullStrip -monkey` *(nhp)* |
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
- **`-S nhp`**: FSL averages the full eddy-corrected (ECC) series into one
  volume, AFNI `3dSkullStrip -monkey` extracts the brain (a preset tuned for
  NHP EPI contrast). The raw mask is used as-is (binarised only — no
  dilation or erosion).

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

```text
process_dti.sh [OPTIONS]

Required:
  -a FILE   AP-encoded DWI series (NIfTI)
  -v FILE   AP bvec file
  -b FILE   AP bval file
  -p FILE   PA-encoded b0 series for topup correction (NIfTI)
  -o DIR    Output directory (created if needed)

Optional:
  -S SPEC   Species / masking: 'human' (dwi2mask) or 'nhp'
            (AFNI 3dSkullStrip -monkey, raw mask)            [default: human]
  -f FILE   FA template (NIfTI) → also warp DTI maps to template via ANTs SyN
  -r FLOAT  Total EPI readout time (s)   [default: 0.06 or JSON sidecar]
  -s INT    DTI shell for tensor fitting [default: 1000]
  -c FILE   Pre-made brain mask (NIfTI) → resume from bias correction onward
  -k        Keep intermediate files in <outdir>/tmp/
  -h        Help
```

> **Invocation forms.** The container's `run` action defaults to `process_dti.sh`, so
> `apptainer run <sif> -a ... -S nhp` (arguments only, no script path) runs the pipeline.
> The examples below use the explicit `apptainer exec <sif> /app/process_dti.sh ...` form.
> Do **not** combine them as `apptainer run <sif> /app/process_dti.sh ...` — the script
> path would be passed to the script as an argument and option parsing would fail.

### Human data

```bash
apptainer exec \
  -B /data/study:/data -B /data/out:/out \
  alternative_dti_processing.sif \
  /app/process_dti.sh \
      -a /data/sub-01_AP.nii.gz \
      -v /data/sub-01_AP.bvec \
      -b /data/sub-01_AP.bval \
      -p /data/sub-01_PA.nii.gz \
      -o /out/sub-01
```

### Monkey / NHP data (AFNI `-monkey` masking)

```bash
apptainer exec \
  -B /data/nhp_study:/data -B /data/out:/out \
  alternative_dti_processing.sif \
  /app/process_dti.sh \
      -a /data/monkey01_AP.nii.gz \
      -v /data/monkey01_AP.bvec \
      -b /data/monkey01_AP.bval \
      -p /data/monkey01_PA.nii.gz \
      -o /out/monkey01 \
      -S nhp
```

### Optional: warp DTI maps to an FA template

Add `-f /data/FA_template.nii.gz` to either invocation. After fitting the
native maps, the pipeline registers native FA → template with ANTs SyN and
applies the composite warp to all four metrics, producing
`w<base>_ECC_{FA,MD,AD,RD}.nii.gz`.

### Resuming with a hand-edited mask

If an auto-generated mask needs manual correction, edit
`<outdir>/<base>_mask.nii.gz` (or supply your own in the same space/resolution
as `<base>_ECC.nii.gz`) and rerun with `-c`:

```bash
/app/process_dti.sh -a ... -v ... -b ... -p ... -o /out/monkey01 \
    -c /out/monkey01/monkey01_AP_mask_edited.nii.gz
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

## Notes

- **Multi-shell data** (b = 0, 500, 1000, 2000): the tensor is fitted on
  b=0 + b=1000 by default (`-s 2000` to use the high shell). The full
  preprocessed series is preserved for advanced models (NODDI, SMT, MSMT-CSD).
- **Phase-encode direction**: STEP 5 assumes AP = `j-` (`-pe_dir AP`). Verify
  with `mrinfo` / your protocol; flip if your AP is actually posterior→anterior.
- **Readout time**: auto-read from the AP JSON sidecar (`TotalReadoutTime`) when
  present; override with `-r`.

## Documentation

- [`methods_nhp.md`](methods_nhp.md) — methods-section write-up of the NHP pipeline.
- [`NHP_DTI_Pipeline_Guide.md`](NHP_DTI_Pipeline_Guide.md) — step-by-step student guide.
