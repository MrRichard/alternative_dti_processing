#!/bin/bash
#
# process_dti.sh — DTI Processing Pipeline for Multi-shell AP/PA Acquisition
#
# Designed to run inside the MRtrix3 Singularity container.
# Requires tools: MRtrix3 (mrconvert, dwidenoise, mrdegibbs, dwifslpreproc,
#                           dwibiascorrect, dwi2mask, dwi2tensor, tensor2metric,
#                           dwiextract, mrmath, mrcat)
#                 FSL (topup/eddy — available in container via dwifslpreproc dependency)
#
# ------------------------------------------------------------------------------------------------------------------------------------------
# INPUT DATA NOTES
# ------------------------------------------------------------------------------------------------------------------------------------------
# This pipeline handles multi-shell DWI with b-values: 0, 500, 1000, 2000.
#   b=0    — unweighted reference volumes
#   b=500  — low-weighting shell (not used for tensor fitting)
#   b=1000 — standard DTI shell (used for FA/MD/AD/RD by default)
#   b=2000 — high-weighting shell (suitable for HARDI/NODDI; not used here)
#
# Tensor fitting uses only b=0 + b=1000 (standard clinical DTI approach).
# The full multi-shell preprocessed data is saved for advanced modeling.
# ------------------------------------------------------------------------------------------------------------------------------------------
#
# Usage (inside container):
#   ./process_dti.sh \
#       -a AP_DWI.nii.gz     \
#       -v AP.bvec           \
#       -b AP.bval           \
#       -p PA_DWI.nii.gz     \
#       -o output_dir

set -euo pipefail

# =============================================================================
# FUNCTIONS
# =============================================================================

usage() {
    cat <<EOF

DTI Processing Pipeline — Multi-shell AP/PA DWI Data

Usage: $(basename "$0") [OPTIONS]

Required:
  -a FILE   AP-encoded DWI series (NIfTI)
  -v FILE   AP bvec file
  -b FILE   AP bval file
  -p FILE   PA-encoded b0 series for topup correction (NIfTI)
  -o DIR    Output directory (created if it does not exist)

Optional:
  -r FLOAT  Total EPI readout time in seconds [default: 0.06]
            If not specified, the value is read from the AP DWI JSON sidecar
            (TotalReadoutTime field) when present, otherwise 0.06 is used.
            Derivation: EffectiveEchoSpacing_s * (PhaseEncodeLines - 1)
            Check your DICOM header or scanner protocol sheet.
  -s INT    B-value shell to use for DTI tensor fitting [default: 1000]
            Use -s 2000 to fit at the high b-value shell instead.
  -c FILE   Path to a pre-made brain mask (NIfTI) in the same space and
            resolution as the ECC NIfTI output. Use this to resume from
            bias correction onward, skipping topup/eddy/auto-masking.
            Mutually exclusive with -m.
  -m        Manual skull-strip mode. The pipeline runs through eddy-current
            correction (step 5), saves the ECC NIfTI+bvec+bval, generates
            a ready-to-run AFNI skull-stripping script, then stops so you
            can create a brain mask manually on a cluster with AFNI installed.
            Mutually exclusive with -c.
  -k        Keep intermediate files in <outdir>/tmp/ [default: delete]
  -h        Show this help message

Outputs (in <outdir>/):
  <apbase>_ECC.{nii.gz,bvec,bval}     — preprocessed full multi-shell DWI
  <apbase>_mask.nii.gz
  <apbase>_{FA,MD,AD,RD}.nii.gz       — DTI maps in native DWI space
  <apbase>_make_mask.sh              — (with -m only) AFNI skull-strip script
  pipeline.log                         — full processing log

NHP workflow (non-human primate data with failed auto-masking):
  Step 1 — Run with -m to generate ECC output and the AFNI mask script:
    ./process_dti.sh -a ... -v ... -b ... -p ... -o OUTDIR -m

  Step 2 — Run the generated AFNI script on a cluster with AFNI installed:
    bash ${OUTDIR}/${AP_BASE}_make_mask.sh
    # This produces: ${OUTDIR}/${AP_BASE}_brain_mask.nii.gz

  Step 3 — Resume with -c pointing to the mask:
    ./process_dti.sh -a ... -v ... -b ... -p ... -o OUTDIR \
                    -c ${OUTDIR}/${AP_BASE}_brain_mask.nii.gz
    (omit -m; all other flags identical to step 1)

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
KEEP_INTERMEDIATES=false
MANUAL_SKULL_STRIP=false
ECC_MASK=""

while getopts "a:v:b:p:o:r:s:c:mkth" opt; do
    case "$opt" in
        a) AP_DWI="$OPTARG" ;;
        v) AP_BVEC="$OPTARG" ;;
        b) AP_BVAL="$OPTARG" ;;
        p) PA_DWI="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        r) READOUT_TIME="$OPTARG"; READOUT_TIME_MANUAL=true ;;
        s) DTI_SHELL="$OPTARG" ;;
        c) ECC_MASK="$OPTARG" ;;
        m) MANUAL_SKULL_STRIP=true ;;
        k) KEEP_INTERMEDIATES=true ;;
        h) usage ;;
        *) die "Unknown option. Use -h for help." ;;
    esac
done

[[ -z "${AP_DWI:-}"  ]] && die "Missing required argument: -a AP_DWI"
[[ -z "${AP_BVEC:-}" ]] && die "Missing required argument: -v AP_BVEC"
[[ -z "${AP_BVAL:-}" ]] && die "Missing required argument: -b AP_BVAL"
[[ -z "${PA_DWI:-}"  ]] && die "Missing required argument: -p PA_DWI"
[[ -z "${OUTDIR:-}"  ]] && die "Missing required argument: -o OUTDIR"

[[ -f "$AP_DWI"  ]] || die "AP DWI not found: $AP_DWI"
[[ -f "$AP_BVEC" ]] || die "AP bvec not found: $AP_BVEC"
[[ -f "$AP_BVAL" ]] || die "AP bval not found: $AP_BVAL"
[[ -f "$PA_DWI"  ]] || die "PA DWI not found: $PA_DWI"

# Validate mask-related arguments
if [[ -n "$ECC_MASK" && ! -f "$ECC_MASK" ]]; then
    die "Pre-made brain mask not found: $ECC_MASK"
fi

# Mutually exclusive: -m and -c
if [[ "$MANUAL_SKULL_STRIP" = true && -n "$ECC_MASK" ]]; then
    die "Options -m and -c are mutually exclusive. Use -m to stop after ECC and create a mask, or -c to supply an existing mask — not both."
fi

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

# =============================================================================
# RESUME DETECTION — manual mask workflow
# =============================================================================
RESUME_FROM_MASK=false
if [[ -n "$ECC_MASK" ]]; then
    ECC_NII="$OUTDIR/${AP_BASE}_ECC.nii.gz"
    ECC_BVEC="$OUTDIR/${AP_BASE}_ECC.bvec"
    ECC_BVAL="$OUTDIR/${AP_BASE}_ECC.bval"

    if [[ ! -f "$ECC_NII" ]]; then
        die "Cannot resume: ${ECC_NII} not found. Run with -m first to generate the ECC output, then create a brain mask and rerun with -c."
    fi
    [[ -f "$ECC_BVEC" ]] || die "Cannot resume: ${ECC_BVEC} not found. Run with -m first."
    [[ -f "$ECC_BVAL" ]] || die "Cannot resume: ${ECC_BVAL} not found. Run with -m first."

    RESUME_FROM_MASK=true
fi

# =============================================================================
# LOG HEADER
# =============================================================================
log "============================================================"
log " DTI Processing Pipeline — Multi-shell AP/PA"
log "============================================================"
log " AP DWI               : $AP_DWI"
log " AP bvec              : $AP_BVEC"
log " AP bval              : $AP_BVAL"
log " PA DWI (topup)       : $PA_DWI"
log " Output directory     : $OUTDIR"
log " EPI readout time     : ${READOUT_TIME}s  [source: ${READOUT_TIME_SOURCE}]"
log " DTI shell            : b=${DTI_SHELL}"
log " Keep intermediates   : $KEEP_INTERMEDIATES"
if [[ "$RESUME_FROM_MASK" = true ]]; then
    log " Resume from mask     : YES (${ECC_MASK})"
elif [[ "$MANUAL_SKULL_STRIP" = true ]]; then
    log " Manual skull-strip   : YES (pipeline will stop after ECC)"
fi
log "------------------------------------------------------------"

# Verify required tools are available
log "Checking required tools..."
for tool in mrconvert mrinfo dwidenoise mrdegibbs dwifslpreproc \
            dwibiascorrect dwi2mask dwi2tensor tensor2metric dwiextract \
            mrmath mrcat; do
    command -v "$tool" >/dev/null 2>&1 \
        || die "Required tool not found in PATH: $tool"
done
log "All required tools found."

# =============================================================================
# BRANCH: RESUME FROM MASK vs. NORMAL PIPELINE
# =============================================================================
if [[ "$RESUME_FROM_MASK" = true ]]; then

    # -----------------------------------------------------------------------
    # RESUME PATH: re-import existing ECC output + user-supplied mask
    # -----------------------------------------------------------------------
    log "------------------------------------------------------------"
    log ">> RESUME PATH: importing existing ECC outputs and user-supplied mask"

    ECC_NII="$OUTDIR/${AP_BASE}_ECC.nii.gz"
    ECC_BVEC="$OUTDIR/${AP_BASE}_ECC.bvec"
    ECC_BVAL="$OUTDIR/${AP_BASE}_ECC.bval"

    # Re-import ECC NIfTI+bvec+bval into MIF
    mrconvert "$ECC_NII" \
        -fslgrad "$ECC_BVEC" "$ECC_BVAL" \
        "$TMPDIR/ap_preproc.mif" -force

    log " Imported ECC NIfTI dimensions : $(mrinfo -size "$TMPDIR/ap_preproc.mif")"

    # Import user-supplied brain mask into MIF
    mrconvert "$ECC_MASK" "$TMPDIR/brain_mask.mif" -force

    ECC_VOX=$(mrinfo -size "$TMPDIR/ap_preproc.mif" | awk '{print $1, $2, $3}')
    MASK_VOX=$(mrinfo -size "$TMPDIR/brain_mask.mif" | awk '{print $1, $2, $3}')

    log " ECC voxel dimensions (LxLyLz) : $ECC_VOX"
    log " Mask voxel dimensions (LxLyLz) : $MASK_VOX"

    if [[ "$MASK_VOX" != "$ECC_VOX" ]]; then
        die "Mask dimensions ($MASK_VOX) do not match ECC image dimensions ($ECC_VOX). " \
            "The mask must be in the same space and resolution as the ECC NIfTI."
    fi

    log " Mask dimensions verified — matches ECC image."
    log " Saved: ${AP_BASE}_mask.nii.gz (copy of user-supplied mask)"

    # Copy user mask to output as well
    mrconvert "$TMPDIR/brain_mask.mif" "$OUTDIR/${AP_BASE}_mask.nii.gz" -force

else

    # =============================================================================
    # NORMAL PATH: steps 1-6 (raw → preprocessed + auto-mask)
    # =============================================================================

    # -----------------------------------------------------------------------
    # STEP 1: Convert inputs to MIF format
    # -----------------------------------------------------------------------
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

    # -----------------------------------------------------------------------
    # STEP 2: Denoise — Marchenko-Pastur PCA
    # -----------------------------------------------------------------------
    log "------------------------------------------------------------"
    log "STEP 2: Denoising (MP-PCA)"
    log " NOTE: Denoising must be applied before any other spatial operations."

    dwidenoise \
        "$TMPDIR/ap_raw.mif" \
        "$TMPDIR/ap_denoised.mif" \
        -noise "$TMPDIR/noise_map.mif" \
        -force

    # -----------------------------------------------------------------------
    # STEP 3: Gibbs ringing correction
    # -----------------------------------------------------------------------
    log "------------------------------------------------------------"
    log "STEP 3: Gibbs ringing correction"
    log " NOTE: Applied after denoising, before geometric corrections."

    mrdegibbs \
        "$TMPDIR/ap_denoised.mif" \
        "$TMPDIR/ap_degibbs.mif" \
        -force

    # -----------------------------------------------------------------------
    # STEP 4: Build SE-EPI b0 pair for topup
    # -----------------------------------------------------------------------
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

    # -----------------------------------------------------------------------
    # STEP 5: Distortion + eddy current correction
    # -----------------------------------------------------------------------
    log "------------------------------------------------------------"
    log "STEP 5: Distortion and eddy current correction"
    log " Uses FSL topup (susceptibility) + eddy (eddy currents, motion)"
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

    # -----------------------------------------------------------------------
    # EARLY EXIT: manual skull-strip mode — stop here, let user create mask
    # -----------------------------------------------------------------------
    if [[ "$MANUAL_SKULL_STRIP" = true ]]; then
        MASK_SCRIPT="$OUTDIR/${AP_BASE}_make_mask.sh"

        # Generate the AFNI skull-stripping script for the user to run on a cluster
        cat > "$MASK_SCRIPT" << 'AFNISCRIPT'
#!/bin/bash
# =============================================================================
# AFNI NHP Brain Mask Generation Script
# Generated by process_dti.sh — manual skull-strip workflow
# =============================================================================
#
# HOW TO USE:
#   Run this script on a machine/cluster with AFNI installed.
#   It reads the ECC NIfTI from the output directory and produces a brain mask.
#
#   IMPORTANT: This script must be run AFTER running process_dti.sh with -m.
#   Do NOT move or rename the ECC NIfTI — the mask must be in the same
#   space and resolution as the ECC image.
#
# OUTPUT:
#   <outdir>/<apbase>_brain_mask.nii.gz  — binary brain mask (1=inside, 0=outside)
#
# WHAT THIS SCRIPT DOES:
#   1. Computes the mean of the 4D ECC volume (average across all volumes)
#   2. Runs AFNI 3dSkullStrip with -monkey flag (designed for NHP EPI data)
#   3. Erodes the mask by 1 voxel to remove skull/dura edge signals
#   4. Fills interior holes (ventricles, cysts)
#   5. Keeps only the largest connected component (removes eye blobs, neck)
#
# AFNI INSTALLATION (if not already available):
#   On Ubuntu/Debian:  apt-get install afni
#   On CentOS/RHEL:    yum install afni
#   Or download:       https://afni.nimh.nih.gov/pub/dist/tgz/linux_afni.tgz
#
# =============================================================================

set -e

# --- Configuration (edit these if your setup differs) ---
# These are auto-populated by process_dti.sh based on your data.

AFNISCRIPT

        # Append the actual runtime values into the generated script
        cat >> "$MASK_SCRIPT" << EOF
ECC_NII="${OUTDIR}/${AP_BASE}_ECC.nii.gz"
MASK_OUT="${OUTDIR}/${AP_BASE}_brain_mask.nii.gz"
TMP_DIR=$(mktemp -d)

echo "AFNI NHP Mask Generation"
echo "========================"
echo "Input ECC  : \$ECC_NII"
echo "Output mask: \$MASK_OUT"
echo ""

# Verify input exists
[[ -f "\$ECC_NII" ]] || { echo "ERROR: ECC NIfTI not found: \$ECC_NII"; echo "Run process_dti.sh with -m first."; exit 1; }

# Step 1 — Compute 3D mean of the 4D ECC volume
echo "[1/4] Computing mean of ECC volume..."
3dTstat -mean -prefix "\${TMP_DIR}/ecc_mean.nii.gz" "\$ECC_NII" -overwrite

# Step 2 — Run 3dSkullStrip with -monkey flag (designed for NHP EPI)
echo "[2/4] Running AFNI 3dSkullStrip -monkey..."
3dSkullStrip \
    -input  "\${TMP_DIR}/ecc_mean.nii.gz" \
    -prefix "\${TMP_DIR}/mask_raw.nii.gz" \
    -monkey \
    -mask_vol yes \
    -overwrite

# Step 3 — Erode slightly to remove skull/dura edge signal bleed
echo "[3/4] Eroding mask (1 voxel) and filling holes..."
3dmask_tool \
    -input  "\${TMP_DIR}/mask_raw.nii.gz" \
    -prefix "\${TMP_DIR}/mask_eroded.nii.gz" \
    -erode  1 \
    -fill_holes \
    -overwrite

# Step 4 — Keep only the largest connected component (removes eye blobs, neck)
echo "[4/4] Keeping largest connected component..."
3dmask_tool \
    -input  "\${TMP_DIR}/mask_eroded.nii.gz" \
    -prefix "\$MASK_OUT" \
    -infracluster 0.9 \
    -overwrite

# Verify output
if [[ -f "\$MASK_OUT" ]]; then
    echo ""
    echo "SUCCESS: Brain mask saved to:"
    echo "  \$MASK_OUT"
    echo ""
    echo "Next step — resume the DTI pipeline with:"
    echo "  ./process_dti.sh \\"
    echo "      -a $AP_DWI \\"
    echo "      -v $AP_BVEC \\"
    echo "      -b $AP_BVAL \\"
    echo "      -p $PA_DWI \\"
    echo "      -o $OUTDIR \\"
    echo "      -c \$MASK_OUT"
else
    echo "ERROR: Mask was not produced. Check AFNI installation and try again."
    exit 1
fi

rm -rf "\${TMP_DIR}"
EOF

        chmod +x "$MASK_SCRIPT"

        log ""
        log "============================================================"
        log " MANUAL SKULL-STRIP MODE — pipeline stopped after ECC"
        log "============================================================"
        log ""
        log " The ECC NIfTI and gradient files have been saved:"
        log "   ${OUTDIR}/${AP_BASE}_ECC.nii.gz"
        log "   ${OUTDIR}/${AP_BASE}_ECC.bvec"
        log "   ${OUTDIR}/${AP_BASE}_ECC.bval"
        log ""
        log " A ready-to-run AFNI skull-stripping script has been generated:"
        log "   ${MASK_SCRIPT}"
        log ""
        log " NEXT STEPS:"
        log "  1. Copy the generated script to a machine/cluster with AFNI installed."
        log "  2. Run:  bash ${AP_BASE}_make_mask.sh"
        log "     (This produces: ${OUTDIR}/${AP_BASE}_brain_mask.nii.gz)"
        log "  3. Copy the mask back to this machine if needed."
        log "  4. Rerun this script with identical -a/-v/-b/-p/-o arguments"
        log "     plus -c pointing to the mask (omit -m):"
        log ""
        log "   ./process_dti.sh \\"
        log "       -a $AP_DWI \\"
        log "       -v $AP_BVEC \\"
        log "       -b $AP_BVAL \\"
        log "       -p $PA_DWI \\"
        log "       -o $OUTDIR \\"
        log "       -c ${OUTDIR}/${AP_BASE}_brain_mask.nii.gz"
        log ""
        log "============================================================"
        exit 0
    fi

    # -----------------------------------------------------------------------
    # STEP 6: Create brain mask — MRtrix3 dwi2mask (standard path)
    # -----------------------------------------------------------------------
    log "------------------------------------------------------------"
    log "STEP 6: Creating brain mask from preprocessed DWI"

    dwi2mask "$TMPDIR/ap_preproc.mif" "$TMPDIR/brain_mask.mif" -force
    mrconvert "$TMPDIR/brain_mask.mif" "$OUTDIR/${AP_BASE}_mask.nii.gz" -force

    log " Saved: ${AP_BASE}_mask.nii.gz"

fi
# =============================================================================
# END OF RESUME / NORMAL BRANCH
# Steps 7–10 operate on $TMPDIR/ap_preproc.mif and $TMPDIR/brain_mask.mif
# regardless of which branch populated them.
# =============================================================================

# -----------------------------------------------------------------------
# STEP 7: Bias field correction
# -----------------------------------------------------------------------
log "------------------------------------------------------------"
log "STEP 7: Bias field correction (FSL FAST)"
log " If ANTs is available in your container, consider dwibiascorrect ants"
log " for superior N4 correction."

dwibiascorrect fsl \
    "$TMPDIR/ap_preproc.mif" \
    "$TMPDIR/ap_biascorr.mif" \
    -mask "$TMPDIR/brain_mask.mif" \
    -force

# -----------------------------------------------------------------------
# STEP 8: Extract DTI shell (b=0 + b=<DTI_SHELL>)
# -----------------------------------------------------------------------
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

# -----------------------------------------------------------------------
# STEP 9: Fit diffusion tensor
# -----------------------------------------------------------------------
log "------------------------------------------------------------"
log "STEP 9: Fitting diffusion tensor (b=${DTI_SHELL} shell)"

dwi2tensor \
    "$TMPDIR/dwi_dti_shell.mif" \
    "$TMPDIR/dt.mif" \
    -mask "$TMPDIR/brain_mask.mif" \
    -force

# -----------------------------------------------------------------------
# STEP 10: Extract DTI metrics in native DWI space
# -----------------------------------------------------------------------
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
log "   ${AP_BASE}_mask.nii.gz"
log ""
log " DTI metrics — Native DWI space:"
log "   ${AP_BASE}_FA.nii.gz     Fractional Anisotropy"
log "   ${AP_BASE}_MD.nii.gz     Mean Diffusivity"
log "   ${AP_BASE}_AD.nii.gz     Axial Diffusivity"
log "   ${AP_BASE}_RD.nii.gz     Radial Diffusivity"
log ""
log " Multi-shell note:"
log "   Data contains b = 0, 500, 1000, 2000."
log "   Tensor fitted using b=0 + b=${DTI_SHELL} only (standard DTI)."
log "   For advanced models (NODDI, SMT, MSMT-CSD), use dwi_preproc.nii.gz"
log "   with all shells."