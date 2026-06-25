# Methods — NHP DTI Processing Pipeline

This document is a methods-section write-up of the non-human primate (NHP)
diffusion tensor imaging (DTI) preprocessing and tensor-estimation pipeline
implemented in `process_dti.sh` (invoked with `-S nhp`). It is intended to be
adapted for the methods section of a manuscript. Tool-specific parameters match
the script exactly; adjust the prose to your acquisition.

## Summary

Diffusion-weighted images were processed with a containerized pipeline built on
MRtrix3 (Tournier et al., 2019), FSL (Smith et al., 2004), AFNI (Cox, 1996), and
ANTs (Avants et al., 2008). All steps were executed inside a single
Apptainer/Singularity image to fix software versions and ensure reproducibility.

## Acquisition assumptions

Multi-shell diffusion data were acquired with anterior–posterior (AP) phase
encoding at b-values of 0, 500, 1000, and 2000 s/mm², together with a short
reverse phase-encoded (posterior–anterior, PA) b = 0 series used for
susceptibility-distortion estimation. The effective total readout time was read
from the BIDS JSON sidecar (`TotalReadoutTime`) when available, or specified
manually.

## Preprocessing

1. **Denoising.** Thermal noise was removed from the raw AP series using
   Marchenko–Pastur PCA denoising (`dwidenoise`; Veraart et al., 2016), applied
   before any other spatial operation.

2. **Gibbs-ringing removal.** Gibbs ringing artifacts were suppressed using the
   local subvoxel-shifts method (`mrdegibbs`; Kellner et al., 2016).

3. **Susceptibility and eddy-current correction.** An SE-EPI b0 pair was
   assembled from the mean AP b0 (after denoising and Gibbs correction) and the
   mean PA b0 volume. Susceptibility-induced distortions, eddy-current
   distortions, and inter-volume head motion were corrected jointly using FSL
   `topup` and `eddy` (Andersson & Sotiropoulos, 2016) orchestrated by the
   MRtrix3 wrapper `dwifslpreproc` (`-rpe_pair`, `-pe_dir AP`,
   `--slm=linear --data_is_shelled`). Diffusion gradient directions were rotated
   to follow the motion correction. The corrected multi-shell series (ECC) was
   retained with its rotated gradient table.

## Brain extraction (NHP-specific)

Reliable brain masking of monkey EPI data could not be achieved with MRtrix3
`dwi2mask` or FSL `bet`; both produced over- or under-inclusive masks on this
data. Instead, a mask was derived directly from the corrected diffusion data:

1. The eddy-corrected b=0 volumes were extracted and averaged into a single
   3D volume (`dwiextract -bzero` + `mrmath mean`), yielding a high-SNR image
   with clear brain/background contrast optimized for EPI.
2. The mean b=0 volume was skull-stripped with AFNI `3dSkullStrip` using the
   `-monkey` preset, which is tuned for non-human primate EPI contrast
   (`-mask_vol`).
3. The resulting mask was binarised and eroded by one voxel (`fslmaths -bin -ero`)
   to remove a residual non-brain halo (dura/CSF/skull edge) that the `-monkey`
   preset tends to leave. Surface expansion was controlled with `-blur_fwhm 2`
   (Gaussian blur stabilises the boundary on noisy EPI data) and `-no_touchup`
   (prevents reclaiming non-brain edge voxels at the end of the strip).

This single brain mask was used for all subsequent steps (bias-field
estimation, tensor fitting, and metric extraction). Masks were visually
inspected; where necessary a corrected mask was supplied and the pipeline
resumed from bias correction onward.

## Tensor estimation

4. **Bias-field correction.** Low-frequency intensity inhomogeneity was
   corrected within the brain mask using `dwibiascorrect fsl`
   (FSL FAST; Zhang et al., 2001).

5. **Shell selection and tensor fit.** The b = 0 and b = 1000 s/mm² volumes were
   extracted (`dwiextract`) and the diffusion tensor was fitted by iteratively
   reweighted linear least squares within the brain mask (`dwi2tensor`).

6. **Scalar maps.** Fractional anisotropy (FA), mean diffusivity (MD), axial
   diffusivity (AD), and radial diffusivity (RD) were computed from the tensor
   (`tensor2metric`) in native diffusion space.

## Optional template normalization

When a target FA template was provided, native FA maps were registered to the
template with ANTs symmetric normalization (SyN; `antsRegistrationSyN.sh -t s`),
and the resulting composite transform (affine + warp) was applied to all four
scalar maps (`antsApplyTransforms`, linear interpolation) to bring them into
template space for group-level analysis.

## Software versions

| Software | Version | Role |
|----------|---------|------|
| MRtrix3 | `mrtrix3/mrtrix3:latest` base image | conversion, denoising, degibbs, tensor fit, metrics |
| FSL | bundled in MRtrix3 base image | topup, eddy, FAST, `fslmaths` |
| AFNI | `linux_openmp_64` | `3dSkullStrip -monkey` brain extraction |
| ANTs | 2.6.5 | SyN registration, transform application |

The container is built on Debian 12 (bookworm), amd64. AFNI is installed from
the official precompiled `linux_openmp_64` binaries; ANTs is the precompiled
Ubuntu 22.04 build (compatible with Debian 12 glibc 2.36).

## References

- Andersson JLR, Sotiropoulos SN (2016). *NeuroImage* 125:1063–1078.
- Avants BB, Epstein CL, Grossman M, Gee JC (2008). *Med Image Anal* 12:26–41.
- Cox RW (1996). *Comput Biomed Res* 29:162–173.
- Kellner E, Dhital B, Kiselev VG, Reisert M (2016). *Magn Reson Med* 76:1574–1581.
- Smith SM, et al. (2004). *NeuroImage* 23:S208–S219.
- Tournier J-D, et al. (2019). *NeuroImage* 202:116137.
- Veraart J, et al. (2016). *NeuroImage* 142:394–406.
- Zhang Y, Brady M, Smith S (2001). *IEEE Trans Med Imaging* 20:45–57.
