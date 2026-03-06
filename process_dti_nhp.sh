#!/bin/bash
#
# process_dti_nhp.sh — DTI Processing Pipeline for NHP Multi-shell AP/PA Acquisition
#
# Adapted from process_dti.sh for non-human primate data.
# Key differences from the human pipeline:
#   - Step 4.5: Registers T1w (moving) -> AP b0 (fixed, affine), then brings
#               the T1w brain mask into DWI space for explicit eddy masking.
#               This avoids dwifslpreproc auto-generating a mask.
#   - Step 11: Adds ANTs SyN normalization to warp native DTI maps to a target
#              FA template for group-level analysis.
#
# Requires tools: MRtrix3 (mrconvert, mrinfo, dwidenoise, mrdegibbs,
#                           dwifslpreproc, dwibiascorrect, dwi2tensor,
#                           tensor2metric, dwiextract, mrmath, mrcat)
#                 FSL (bet, fslmaths -- skull stripping and mask binarisation)
#                 ANTs (antsRegistrationSyN.sh, antsApplyTransforms -- required)
#
# INPUT DATA NOTES
# --------------------------------------------------------------------------
# This pipeline handles multi-shell DWI with b-values: 0, 500, 1000, 2000.
#   b=0    — unweighted reference volumes
#   b=500  — low-weighting shell (not used for tensor fitting)
#   b=1000 — standard DTI shell (used for FA/MD/AD/RD by default)
#   b=2000 — high-weighting shell (suitable for HARDI/NODDI; not used here)
#
# Tensor fitting uses only b=0 + b=1000 (standard clinical DTI approach).
# The full multi-shell preprocessed data is saved for advanced modeling.
# --------------------------------------------------------------------------
#
# Usage (inside container):
#   ./process_dti_nhp.sh \
#       -a AP_DWI.nii.gz     \
#       -v AP.bvec           \
#       -b AP.bval           \
#       -p PA_DWI.nii.gz     \
#       -o output_dir        \
#       -t T1w.nii.gz        \
#       -m T1w_brain_mask.nii.gz \
#       -f FA_template.nii.gz
#

set -euo pipefail

# =============================================================================
# FUNCTIONS
# =============================================================================

usage() {
    cat <<EOF

DTI Processing Pipeline — NHP Multi-shell AP/PA DWI Data

Usage: $(basename "$0") [OPTIONS]

Required:
  -a FILE   AP-encoded DWI series (NIfTI)
  -v FILE   AP bvec file
  -b FILE   AP bval file
  -p FILE   PA-encoded b0 series for topup correction (NIfTI)
  -o DIR    Output directory (created if it does not exist)
  -t FILE   Native T1w image (NIfTI) — used for mask coregistration
  -m FILE   Binary brain mask in T1w native space (NIfTI)
  -f FILE   Target FA template for ANTs SyN normalization (NIfTI)

Optional:
  -r FLOAT  Total EPI readout time in seconds [default: 0.06]
            If not specified, the value is read from the AP DWI JSON sidecar
            (TotalReadoutTime field) when present, otherwise 0.06 is used.
            Derivation: EffectiveEchoSpacing_s * (PhaseEncodeLines - 1)
            Check your DICOM header or scanner protocol sheet.
  -s INT    B-value shell to use for DTI tensor fitting [default: 1000]
            Use -s 2000 to fit at the high b-value shell instead.
  -F FLOAT  BET fractional intensity threshold for b0 skull stripping [default: 0.25]
            Smaller values extract larger brain volume. NHP range: 0.20-0.35.
  -G FLOAT  BET vertical gradient in f (-1 to 1) [default: -0.1]
            Negative shifts centre of gravity inferiorly for NHP brain position.
  -k        Keep intermediate files in <outdir>/tmp/ [default: delete]
  -h        Show this help message

Outputs (in <outdir>/):
  <apbase>_ECC.{nii.gz,bvec,bval}     — preprocessed full multi-shell DWI
  <apbase>_mask.nii.gz                 — brain mask in DWI space (BET on mean corrected b0)
  <apbase>_{FA,MD,AD,RD}.nii.gz       — DTI maps in native DWI space
  w<apbase>_{FA,MD,AD,RD}.nii.gz      — DTI maps warped to FA template space
  pipeline.log                         — full processing log

EOF
    exit 0
}

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOGFILE"
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

# Average a 4D image along the volume axis, or just copy if already 3D
# (mrmath mean -axis 3 fails on 3D input)
average_or_copy() {
    local input="$1"
    local output="$2"
    local ndim
    ndim=$(mrinfo -ndim "$input")
    if [ "$ndim" -ge 4 ]; then
        mrmath "$input" mean -axis 3 "$output" -force
    else
        mrconvert "$input" "$output" -force
    fi
}

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================
READOUT_TIME="0.06"
READOUT_TIME_MANUAL=false
DTI_SHELL="1000"
BET_F="0.25"
BET_G="-0.1"
KEEP_INTERMEDIATES=false

while getopts "a:v:b:p:o:t:m:f:r:s:F:G:kh" opt; do
    case "$opt" in
        a) AP_DWI="$OPTARG" ;;
        v) AP_BVEC="$OPTARG" ;;
        b) AP_BVAL="$OPTARG" ;;
        p) PA_DWI="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) T1W="$OPTARG" ;;
        m) T1W_MASK="$OPTARG" ;;
        f) FA_TEMPLATE="$OPTARG" ;;
        r) READOUT_TIME="$OPTARG"; READOUT_TIME_MANUAL=true ;;
        s) DTI_SHELL="$OPTARG" ;;
        F) BET_F="$OPTARG" ;;
        G) BET_G="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=true ;;
        h) usage ;;
        *) die "Unknown option. Use -h for help." ;;
    esac
done

[[ -z "${AP_DWI:-}"      ]] && die "Missing required argument: -a AP_DWI"
[[ -z "${AP_BVEC:-}"     ]] && die "Missing required argument: -v AP_BVEC"
[[ -z "${AP_BVAL:-}"     ]] && die "Missing required argument: -b AP_BVAL"
[[ -z "${PA_DWI:-}"      ]] && die "Missing required argument: -p PA_DWI"
[[ -z "${OUTDIR:-}"      ]] && die "Missing required argument: -o OUTDIR"
[[ -z "${T1W:-}"         ]] && die "Missing required argument: -t T1W"
[[ -z "${T1W_MASK:-}"    ]] && die "Missing required argument: -m T1W_MASK"
[[ -z "${FA_TEMPLATE:-}" ]] && die "Missing required argument: -f FA_TEMPLATE"

[[ -f "$AP_DWI"      ]] || die "AP DWI not found: $AP_DWI"
[[ -f "$AP_BVEC"     ]] || die "AP bvec not found: $AP_BVEC"
[[ -f "$AP_BVAL"     ]] || die "AP bval not found: $AP_BVAL"
[[ -f "$PA_DWI"      ]] || die "PA DWI not found: $PA_DWI"
[[ -f "$T1W"         ]] || die "T1w image not found: $T1W"
[[ -f "$T1W_MASK"    ]] || die "T1w brain mask not found: $T1W_MASK"
[[ -f "$FA_TEMPLATE" ]] || die "FA template not found: $FA_TEMPLATE"

# Derive base filename from AP DWI input (strip path and .nii/.nii.gz extension)
AP_BASE=$(basename "$AP_DWI")
AP_BASE="${AP_BASE%.nii.gz}"
AP_BASE="${AP_BASE%.nii}"

# Auto-detect TotalReadoutTime from JSON sidecar (unless -r was explicitly set)
READOUT_TIME_SOURCE="default (${READOUT_TIME}s)"
if [ "$READOUT_TIME_MANUAL" = false ]; then
    AP_JSON="$(dirname "$AP_DWI")/${AP_BASE}.json"
    if [[ -f "$AP_JSON" ]]; then
        JSON_TRT=$(python3 -c \
            "import json,sys; d=json.load(open(sys.argv[1])); print(d['TotalReadoutTime'])" \
            "$AP_JSON" 2>/dev/null || true)
        if [[ -n "$JSON_TRT" ]]; then
            READOUT_TIME="$JSON_TRT"
            READOUT_TIME_SOURCE="JSON sidecar (${AP_JSON})"
        else
            READOUT_TIME_SOURCE="default (${READOUT_TIME}s) — TotalReadoutTime not found in ${AP_JSON}"
        fi
    else
        READOUT_TIME_SOURCE="default (${READOUT_TIME}s) — no JSON sidecar found"
    fi
else
    READOUT_TIME_SOURCE="manual -r flag (${READOUT_TIME}s)"
fi

# =============================================================================
# SETUP
# =============================================================================
mkdir -p "$OUTDIR"
TMPDIR="$OUTDIR/tmp"
mkdir -p "$TMPDIR"
LOGFILE="$OUTDIR/pipeline.log"

log "============================================================"
log " DTI Processing Pipeline — NHP Multi-shell AP/PA"
log "============================================================"
log " AP DWI               : $AP_DWI"
log " AP bvec              : $AP_BVEC"
log " AP bval              : $AP_BVAL"
log " PA DWI (topup)       : $PA_DWI"
log " Output directory     : $OUTDIR"
log " T1w image            : $T1W"
log " T1w brain mask       : $T1W_MASK"
log " FA template          : $FA_TEMPLATE"
log " EPI readout time     : ${READOUT_TIME}s  [source: ${READOUT_TIME_SOURCE}]"
log " DTI shell            : b=${DTI_SHELL}"
log " Keep intermediates   : $KEEP_INTERMEDIATES"
log "------------------------------------------------------------"

# Verify required tools are available
log "Checking required tools..."
for tool in mrconvert mrinfo dwidenoise mrdegibbs dwifslpreproc \
            dwibiascorrect dwi2tensor tensor2metric dwiextract \
            mrmath mrcat; do
    command -v "$tool" >/dev/null 2>&1 \
        || die "Required tool not found in PATH: $tool"
done

for tool in bet fslmaths; do
    command -v "$tool" >/dev/null 2>&1 \
        || die "Required tool not found in PATH: $tool (FSL)"
done

command -v antsRegistrationSyN.sh >/dev/null 2>&1 \
    || die "Required tool not found in PATH: antsRegistrationSyN.sh (ANTs)"

command -v antsApplyTransforms >/dev/null 2>&1 \
    || die "Required tool not found in PATH: antsApplyTransforms (ANTs)"

log "All required tools found."

# =============================================================================
# STEP 1: Convert inputs to MIF format
# =============================================================================
log "------------------------------------------------------------"
log "STEP 1: Converting to MIF format"

mrconvert "$AP_DWI" \
    -fslgrad "$AP_BVEC" "$AP_BVAL" \
    "$TMPDIR/ap_raw.mif" -force

mrconvert "$PA_DWI" "$TMPDIR/pa_raw.mif" -force

log " AP DWI dimensions : $(mrinfo -size "$TMPDIR/ap_raw.mif")"
log " PA DWI dimensions : $(mrinfo -size "$TMPDIR/pa_raw.mif")"
log " Detected b-value shells: $(mrinfo -shell_bvalues "$TMPDIR/ap_raw.mif")"
log " Volumes per shell: $(mrinfo -shell_sizes "$TMPDIR/ap_raw.mif")"

# =============================================================================
# STEP 2: Denoise — Marchenko-Pastur PCA
# =============================================================================
log "------------------------------------------------------------"
log "STEP 2: Denoising (MP-PCA)"
log " NOTE: Denoising must be applied before any other spatial operations."

dwidenoise \
    "$TMPDIR/ap_raw.mif" \
    "$TMPDIR/ap_denoised.mif" \
    -noise "$TMPDIR/noise_map.mif" \
    -force

# =============================================================================
# STEP 3: Gibbs ringing correction
# =============================================================================
log "------------------------------------------------------------"
log "STEP 3: Gibbs ringing correction"
log " NOTE: Applied after denoising, before geometric corrections."

mrdegibbs \
    "$TMPDIR/ap_denoised.mif" \
    "$TMPDIR/ap_degibbs.mif" \
    -force

# =============================================================================
# STEP 4: Build SE-EPI b0 pair for topup
# =============================================================================
log "------------------------------------------------------------"
log "STEP 4: Building SE-EPI b0 pair for distortion correction"
log " AP b0 volumes: extracted from degibbed AP series"
log " PA b0 volumes: all volumes from PA series (assumed all b=0)"
log " Order in se_epi: [AP b0, PA b0] — required for dwifslpreproc -rpe_pair"

# Mean AP b0 (after denoising + Gibbs)
dwiextract -bzero "$TMPDIR/ap_degibbs.mif" - | \
    mrmath - mean -axis 3 "$TMPDIR/b0_AP.mif" -force

# Mean PA volumes (the short topup series; all expected to be b=0)
average_or_copy "$TMPDIR/pa_raw.mif" "$TMPDIR/b0_PA.mif"

# Concatenate into a 2-volume 4D SE-EPI image: AP then PA
mrcat "$TMPDIR/b0_AP.mif" "$TMPDIR/b0_PA.mif" \
    -axis 3 "$TMPDIR/se_epi_pair.mif" -force

log " SE-EPI pair dimensions: $(mrinfo -size "$TMPDIR/se_epi_pair.mif")"

# =============================================================================
# STEP 4.5: Build explicit eddy mask from native T1w mask
# =============================================================================
log "------------------------------------------------------------"
log "STEP 4.5: Preparing eddy mask from T1w brain mask"
log " Registering T1w -> AP b0 (affine) and resampling T1w mask into DWI space"

# Convert AP b0 to NIfTI reference for ANTs
mrconvert "$TMPDIR/b0_AP.mif" "$TMPDIR/b0_AP.nii.gz" -force

# Fast affine cross-modality registration for eddy masking
antsRegistrationSyN.sh \
    -d 3 \
    -f "$TMPDIR/b0_AP.nii.gz" \
    -m "$T1W" \
    -o "$TMPDIR/t1_to_apb0_" \
    -t a

# Bring T1w brain mask into AP DWI space with nearest-neighbor interpolation
antsApplyTransforms \
    -d 3 \
    -i "$T1W_MASK" \
    -r "$TMPDIR/b0_AP.nii.gz" \
    -t "$TMPDIR/t1_to_apb0_0GenericAffine.mat" \
    -o "$TMPDIR/eddy_mask_from_t1.nii.gz" \
    -n NearestNeighbor

# Enforce strict binary mask and convert to MIF for dwifslpreproc
fslmaths "$TMPDIR/eddy_mask_from_t1.nii.gz" -bin "$TMPDIR/eddy_mask_from_t1_bin.nii.gz"
mrconvert "$TMPDIR/eddy_mask_from_t1_bin.nii.gz" "$TMPDIR/eddy_mask_from_t1.mif" -force

log " Prepared eddy mask: $TMPDIR/eddy_mask_from_t1_bin.nii.gz"

# =============================================================================
# STEP 5: Distortion + eddy current correction
# =============================================================================
log "------------------------------------------------------------"
log "STEP 5: Distortion and eddy current correction"
log " Uses FSL topup (susceptibility) + eddy (eddy currents, motion)"
log " Using explicit eddy mask from T1w (prevents internal auto-mask generation)"
log " This step typically takes 30–90 minutes depending on data size."
log ""
log " IMPORTANT: -pe_dir AP assumes your AP series is j- phase-encoded"
log " (anterior → posterior). Verify with mrinfo or your scanner protocol."
log " If your AP direction is actually PA (posterior → anterior), swap to -pe_dir PA."

dwifslpreproc \
    "$TMPDIR/ap_degibbs.mif" \
    "$TMPDIR/ap_preproc.mif" \
    -rpe_pair \
    -se_epi "$TMPDIR/se_epi_pair.mif" \
    -pe_dir AP \
    -readout_time "$READOUT_TIME" \
    -eddy_mask "$TMPDIR/eddy_mask_from_t1.mif" \
    -eddy_options "--slm=linear --data_is_shelled" \
    -force

# Export full preprocessed multi-shell DWI (NIfTI + FSL gradients)
# Eddy correction rotates b-vectors — these are the corrected gradients.
mrconvert "$TMPDIR/ap_preproc.mif" \
    "$OUTDIR/${AP_BASE}_ECC.nii.gz" \
    -export_grad_fsl \
        "$OUTDIR/${AP_BASE}_ECC.bvec" \
        "$OUTDIR/${AP_BASE}_ECC.bval" \
    -force

log " Saved preprocessed DWI: ${AP_BASE}_ECC.nii.gz / .bvec / .bval"

# =============================================================================
# STEP 6: Create brain mask — FSL BET on mean corrected b0
# =============================================================================
# NHP-specific: dwi2mask and structural coregistration approaches have proven
# unreliable for this dataset. FSL BET applied directly to the mean corrected
# b0 is used instead.
#
# Why BET on the b0 (not the T1w or the full DWI):
#   - The ECC-corrected mean b0 has good SNR and clear brain/background contrast
#   - Avoids cross-modality registration uncertainty entirely
#   - Faster and more predictable than ANTs SyN at this stage
#
# Why a mask is needed before tensor fitting:
#   - dwibiascorrect uses the mask to estimate the smooth bias field;
#     including skull/soft tissue pulls the estimate off for brain voxels
#   - The FA map fed to ANTs in Step 11 will have artifactual high-FA voxels
#     at skull interfaces without a mask, degrading template registration
#
# NHP BET tuning notes:
#   -f FLOAT  Fractional intensity threshold (0->1); smaller = larger brain
#             extract. NHP typically needs 0.20-0.35 (human default is 0.5).
#   -g FLOAT  Vertical gradient in f (-1->1); negative shifts the centre
#             of gravity inferiorly to compensate for NHP brain position.
#             Try -0.05 to -0.15 as a starting point.
#   -R        Robust iterative re-estimation; recommended when default sphere
#             initialisation is poor (common in NHP).
#   -m        Output binary mask file (<out>_mask.nii.gz).
#
# These defaults (-f 0.25 -g -0.1 -R) are a reasonable NHP starting point
# but MUST be inspected and tuned per dataset. Use the -F and -G pipeline
# flags to override without editing this script.
# =============================================================================
log "------------------------------------------------------------"
log "STEP 6: Creating brain mask via FSL BET on mean b0 (NHP)"
log " BET parameters: -f ${BET_F} -g ${BET_G} -R -m"
log " Override defaults with pipeline flags -F (f threshold) and -G (gradient)"

# 6a: Extract mean corrected b0 from preprocessed DWI
log " 6a: Extracting mean corrected b0"
dwiextract -bzero "$TMPDIR/ap_preproc.mif" - | \
    mrmath - mean -axis 3 "$TMPDIR/b0_corrected.mif" -force
mrconvert "$TMPDIR/b0_corrected.mif" "$TMPDIR/b0_corrected.nii.gz" -force

# 6b: Run BET on the mean b0
#     -f and -g are set by pipeline flags -F/-G (defaults: 0.25 / -0.1)
#     -R: robust iterative re-estimation (recommended for NHP)
#     -m: write binary mask as b0_bet_mask.nii.gz
log " 6b: Running BET"
bet \
    "$TMPDIR/b0_corrected.nii.gz" \
    "$TMPDIR/b0_bet" \
    -f "$BET_F" \
    -g "$BET_G" \
    -R \
    -m

# 6c: Binarise mask (BET -m output should already be binary, but enforce it
#     to guard against any float encoding from older FSL versions)
log " 6c: Binarising mask"
fslmaths "$TMPDIR/b0_bet_mask.nii.gz" -bin "$TMPDIR/brain_mask_DWI.nii.gz"

# 6d: Convert mask to MIF for downstream MRtrix steps
mrconvert "$TMPDIR/brain_mask_DWI.nii.gz" "$TMPDIR/brain_mask.mif" -force

mrconvert "$TMPDIR/brain_mask.mif" "$OUTDIR/${AP_BASE}_mask.nii.gz" -force
log " Saved: ${AP_BASE}_mask.nii.gz"
log " QC: inspect $TMPDIR/b0_bet.nii.gz and ${AP_BASE}_mask.nii.gz"
log "     If mask is over/under-inclusive, re-run with adjusted -F and -G flags"


# =============================================================================
# STEP 7: Bias field correction
# =============================================================================
log "------------------------------------------------------------"
log "STEP 7: Bias field correction (FSL FAST)"
log " If ANTs is available in your container, consider dwibiascorrect ants"
log " for superior N4 correction."

dwibiascorrect fsl \
    "$TMPDIR/ap_preproc.mif" \
    "$TMPDIR/ap_biascorr.mif" \
    -mask "$TMPDIR/brain_mask.mif" \
    -force

# =============================================================================
# STEP 8: Extract DTI shell (b=0 + b=<DTI_SHELL>)
# =============================================================================
log "------------------------------------------------------------"
log "STEP 8: Extracting b=0 + b=${DTI_SHELL} shell for tensor fitting"
log ""
log " Multi-shell data summary:"
log "   b=0    — $(mrinfo -shell_sizes "$TMPDIR/ap_biascorr.mif" | awk 'NR==1{print $1}') reference volumes"
log "   b=500  — low-weighting shell (not used in tensor fit)"
log "   b=1000 — standard DTI shell"
log "   b=2000 — high-weighting shell (HARDI/NODDI; not used in tensor fit)"
log " Fitting tensor using: b=0 + b=${DTI_SHELL}"

dwiextract \
    -shells "0,${DTI_SHELL}" \
    "$TMPDIR/ap_biascorr.mif" \
    "$TMPDIR/dwi_dti_shell.mif" \
    -force

log " DTI shell dimensions: $(mrinfo -size "$TMPDIR/dwi_dti_shell.mif")"

# =============================================================================
# STEP 9: Fit diffusion tensor
# =============================================================================
log "------------------------------------------------------------"
log "STEP 9: Fitting diffusion tensor (b=${DTI_SHELL} shell)"

dwi2tensor \
    "$TMPDIR/dwi_dti_shell.mif" \
    "$TMPDIR/dt.mif" \
    -mask "$TMPDIR/brain_mask.mif" \
    -force

# =============================================================================
# STEP 10: Extract DTI metrics in native DWI space
# =============================================================================
log "------------------------------------------------------------"
log "STEP 10: Extracting DTI metrics (native DWI space)"

tensor2metric \
    "$TMPDIR/dt.mif" \
    -fa  "$TMPDIR/fa.mif" \
    -adc "$TMPDIR/md.mif" \
    -ad  "$TMPDIR/ad.mif" \
    -rd  "$TMPDIR/rd.mif" \
    -mask "$TMPDIR/brain_mask.mif" \
    -force

for metric in fa md ad rd; do
    METRIC_UPPER=$(echo "$metric" | tr '[:lower:]' '[:upper:]')
    mrconvert "$TMPDIR/${metric}.mif" "$OUTDIR/${AP_BASE}_${METRIC_UPPER}.nii.gz" -force
    log " Saved: ${AP_BASE}_${METRIC_UPPER}.nii.gz"
done

# =============================================================================
# STEP 11: Normalize DTI maps to template space (ANTs SyN)
# =============================================================================
log "------------------------------------------------------------"
log "STEP 11: Normalizing DTI maps to FA template space (ANTs SyN)"

# 11a: Register native FA → FA template using SyN
log " 11a: Running antsRegistrationSyN.sh (FA → template)"
log " This step may take 30–60 minutes."
antsRegistrationSyN.sh \
    -d 3 \
    -f "$FA_TEMPLATE" \
    -m "$OUTDIR/${AP_BASE}_FA.nii.gz" \
    -o "$TMPDIR/ants_" \
    -t s

log " ANTs registration complete."
log "   Affine:       $TMPDIR/ants_0GenericAffine.mat"
log "   Warp:         $TMPDIR/ants_1Warp.nii.gz"
log "   Inverse warp: $TMPDIR/ants_1InverseWarp.nii.gz"

# 11b: Apply composite warp to all four DTI metrics
log " 11b: Applying composite warp to DTI metrics"
for metric in fa md ad rd; do
    METRIC_UPPER=$(echo "$metric" | tr '[:lower:]' '[:upper:]')
    NATIVE_MAP="$OUTDIR/${AP_BASE}_${METRIC_UPPER}.nii.gz"
    WARPED_MAP="$OUTDIR/w${AP_BASE}_ECC_${METRIC_UPPER}.nii.gz"

    antsApplyTransforms \
        -d 3 \
        -i "$NATIVE_MAP" \
        -r "$FA_TEMPLATE" \
        -t "$TMPDIR/ants_1Warp.nii.gz" \
        -t "$TMPDIR/ants_0GenericAffine.mat" \
        -o "$WARPED_MAP" \
        -n Linear

    log " Saved (template space): w${AP_BASE}_ECC_${METRIC_UPPER}.nii.gz"
done

# =============================================================================
# CLEANUP
# =============================================================================
log "------------------------------------------------------------"
if [ "$KEEP_INTERMEDIATES" = false ]; then
    log "Removing intermediate files (use -k to retain them)"
    rm -rf "$TMPDIR"
else
    log "Intermediate files retained in: $TMPDIR"
fi

# =============================================================================
# SUMMARY
# =============================================================================
log "============================================================"
log " Pipeline complete!"
log ""
log " Outputs in: $OUTDIR"
log ""
log " Preprocessed DWI (full multi-shell, native space):"
log "   ${AP_BASE}_ECC.nii.gz   eddy-corrected DWI (all shells)"
log "   ${AP_BASE}_ECC.bvec     eddy-rotated gradient directions"
log "   ${AP_BASE}_ECC.bval     b-values" 
log "   ${AP_BASE}_mask.nii.gz  brain mask (BET on mean corrected b0)"
log ""
log " DTI metrics — Native DWI space:"
log "   ${AP_BASE}_FA.nii.gz     Fractional Anisotropy"
log "   ${AP_BASE}_MD.nii.gz     Mean Diffusivity"
log "   ${AP_BASE}_AD.nii.gz     Axial Diffusivity"
log "   ${AP_BASE}_RD.nii.gz     Radial Diffusivity"
log ""
log " DTI metrics — Template space (ANTs SyN):"
log "   w${AP_BASE}_ECC_FA.nii.gz"
log "   w${AP_BASE}_ECC_MD.nii.gz"
log "   w${AP_BASE}_ECC_AD.nii.gz"
log "   w${AP_BASE}_ECC_RD.nii.gz"
log ""
log " Multi-shell note:"
log "   Data contains b = 0, 500, 1000, 2000."
log "   Tensor fitted using b=0 + b=${DTI_SHELL} only (standard DTI)."
log "   For advanced models (NODDI, SMT, MSMT-CSD), use ${AP_BASE}_ECC.nii.gz"
log "   with all shells."
log "============================================================"
