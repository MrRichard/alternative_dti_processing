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
# Optional NHP brain extraction (flag -M):
#   afni    — AFNI 3dSkullStrip with -monkey flag; requires AFNI in container
#   hd-bet  — Deep learning HD-BET; requires Python + hd-bet package installed
#   deepbet — DeepBet U-Net for NHP; requires Python + deepbet package installed
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
#

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
            Mutually exclusive with -m and -M.
  -m        Manual skull-strip mode. The pipeline runs through eddy-current
            correction (step 5), saves the ECC NIfTI+bvec+bval, then stops
            so you can create a brain mask manually (e.g. ITK-SNAP, FSLeyes).
            Rerun with -c /path/to/mask.nii.gz to resume from bias correction.
            Mutually exclusive with -c and -M.
  -M METHOD Automated NHP brain extraction — runs immediately after ECC export,
            replacing dwi2mask. Generates the 3D mean of the ECC volume, runs
            the chosen method to produce a brain mask, and continues the pipeline.
            Valid METHOD values:
              afni    — AFNI 3dSkullStrip with -monkey flag (recommended for EPI)
              hd-bet  — HD-BET deep learning brain extraction
              deepbet — DeepBet U-Net trained on NHP data
            Mutually exclusive with -c and -m.
  -k        Keep intermediate files in <outdir>/tmp/ [default: delete]
  -h        Show this help message

Outputs (in <outdir>/):
  <apbase>_ECC.{nii.gz,bvec,bval}     — preprocessed full multi-shell DWI
  <apbase>_mask.nii.gz
  <apbase>_{FA,MD,AD,RD}.nii.gz       — DTI maps in native DWI space
  pipeline.log                         — full processing log

NHP workflow (non-human primate data with failed auto-masking):
  Option A — Automated (no manual mask creation):
    ./process_dti.sh -a ... -v ... -b ... -p ... -o OUTDIR -M afni
    (substitute afni with hd-bet or deepbet if preferred)
  Option B — Manual mask creation:
    1. Run:  ./process_dti.sh -a ... -v ... -b ... -p ... -o OUTDIR -m
    2. Create a brain mask for the ECC NIfTI in ITK-SNAP/FSLeyes/3D Slicer
       (draw/paint over the brain, save as NIfTI in the same space/resolution).
    3. Resume: ./process_dti.sh -a ... -v ... -b ... -p ... -o OUTDIR \\
                    -c /path/to/your/mask.nii.gz
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
NHP_MASK_METHOD=""

while getopts "a:v:b:p:o:r:s:c:mM:kh" opt; do
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
        M) NHP_MASK_METHOD="$OPTARG" ;;
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

# Validate NHP mask method
if [[ -n "$NHP_MASK_METHOD" ]]; then
    case "$NHP_MASK_METHOD" in
        afni|hd-bet|deepbet) ;;
        *) die "Unknown -M method: '$NHP_MASK_METHOD'. Valid options: afni, hd-bet, deepbet" ;;
    esac
fi

# Validate mask-related arguments
if [[ -n "$ECC_MASK" && ! -f "$ECC_MASK" ]]; then
    die "Pre-made brain mask not found: $ECC_MASK"
fi

# Mutually exclusive: -m, -c, -M
MASK_MODES=0
[[ "$MANUAL_SKULL_STRIP" = true ]] && ((MASK_MODES++))
[[ -n "$ECC_MASK" ]] && ((MASK_MODES++))
[[ -n "$NHP_MASK_METHOD" ]] && ((MASK_MODES++))

if [[ "$MASK_MODES" -gt 1 ]]; then
    die "Options -m, -c, and -M are mutually exclusive. Choose only one masking approach."
fi

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
elif [[ -n "$NHP_MASK_METHOD" ]]; then
    log " NHP auto-mask method : $NHP_MASK_METHOD"
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

# Verify NHP masking tool is available
if [[ -n "$NHP_MASK_METHOD" ]]; then
    case "$NHP_MASK_METHOD" in
        afni)
            command -v 3dSkullStrip >/dev/null 2>&1 \
                || die "AFNI tool '3dSkullStrip' not found in PATH. " \
                       "Either install AFNI or use a different -M method."
            command -v 3dmask_tool >/dev/null 2>&1 \
                || die "AFNI tool '3dmask_tool' not found in PATH."
            log "NHP mask tool (AFNI 3dSkullStrip -monkey) — available."
            ;;
        hd-bet)
            command -v hd-bet >/dev/null 2>&1 \
                || die "HD-BET tool 'hd-bet' not found in PATH. " \
                       "Install with: pip install hd-bet"
            log "NHP mask tool (HD-BET) — available."
            ;;
        deepbet)
            command -v deepbet >/dev/null 2>&1 \
                || die "DeepBet tool 'deepbet' not found in PATH. " \
                       "Install from: https://github.com/HumanBrainED/NHP-BrainExtraction"
            log "NHP mask tool (DeepBet) — available."
            ;;
    esac
fi

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
        log ""
        log "============================================================"
        log " MANUAL SKULL-STRRIP MODE — pipeline stopped after ECC"
        log "============================================================"
        log ""
        log " The ECC NIfTI and gradient files have been saved:"
        log "   ${OUTDIR}/${AP_BASE}_ECC.nii.gz"
        log "   ${OUTDIR}/${AP_BASE}_ECC.bvec"
        log "   ${OUTDIR}/${AP_BASE}_ECC.bval"
        log ""
        log " NEXT STEPS:"
        log "  1. Open ${AP_BASE}_ECC.nii.gz in ITK-SNAP, FSLeyes, or 3D Slicer."
        log "  2. Create a brain mask (draw/paint over the brain; erase skull, "
        log "     dura, and any signal voids outside the brain)."
        log "  3. Save the mask as a NIfTI file in the same space and resolution"
        log "     as the ECC NIfTI — do NOT resample or transform it."
        log "  4. Rerun this script with identical -a/-v/-b/-p/-o arguments"
        log "     plus -c /path/to/your/mask.nii.gz (omit -m):"
        log ""
        log "   ./process_dti.sh \\"
        log "       -a $AP_DWI \\"
        log "       -v $AP_BVEC \\"
        log "       -b $AP_BVAL \\"
        log "       -p $PA_DWI \\"
        log "       -o $OUTDIR \\"
        log "       -c /path/to/your/mask.nii.gz"
        log ""
        log " Intermediates in $TMPDIR are NOT cleaned up on this exit path"
        log " (nice-to-have for debugging; not required for the resume path)."
        log "============================================================"
        # Exit cleanly without cleaning up $TMPDIR
        exit 0
    fi

    # -----------------------------------------------------------------------
    # STEP 6: Create brain mask — NHP automated method OR MRtrix dwi2mask
    # -----------------------------------------------------------------------
    log "------------------------------------------------------------"

    if [[ -n "$NHP_MASK_METHOD" ]]; then

        # =================================================================
        # NHP AUTO-MASK: generate mean → run chosen method → save mask
        # =================================================================
        log "STEP 6: NHP automated brain extraction (method: $NHP_MASK_METHOD)"

        # Compute 3D mean of the ECC 4D volume for skull-stripping
        mrmath "$TMPDIR/ap_preproc.mif" mean -axis 3 "$TMPDIR/ecc_mean_3d.nii.gz" -force
        log " Computed ECC mean volume for skull-stripping input"

        NHP_MASK_NII="$TMPDIR/nhp_brain_mask_raw.nii.gz"

        case "$NHP_MASK_METHOD" in
            afni)
                log " Running AFNI 3dSkullStrip -monkey..."
                3dSkullStrip \
                    -input "$TMPDIR/ecc_mean_3d.nii.gz" \
                    -prefix "$NHP_MASK_NII" \
                    -monkey \
                    -mask_vol yes \
                    -overwrite

                # Erode slightly and clean up isolated blobs (eyes, neck)
                log " Post-processing mask (erode + connected component cleanup)..."
                3dmask_tool \
                    -input "$NHP_MASK_NII" \
                    -prefix "$TMPDIR/nhp_brain_mask_eroded.nii.gz" \
                    -erode 1 \
                    -fill_holes \
                    -infracluster 0.9 \
                    -overwrite

                # Rename eroded mask as the final mask
                mv "$TMPDIR/nhp_brain_mask_eroded.nii.gz" "$NHP_MASK_NII"
                ;;

            hd-bet)
                log " Running HD-BET..."
                hd-bet \
                    -i "$TMPDIR/ecc_mean_3d.nii.gz" \
                    -o "$NHP_MASK_NII" \
                    -device cpu \
                    -mode fast \
                    -tta 0 \
                    -overwrite

                # HD-BET outputs a soft-max probability image; threshold to binary
                # Use the 50th percentile of the brain region as a conservative threshold
                log " Thresholding HD-BET output to binary mask..."
                local THR
                THR=$(3dmaskave -quiet -mask "$NHP_MASK_NII" "$NHP_MASK_NII" 2>/dev/null \
                      | awk '{print $1 * 0.5}')
                3dcalc \
                    -a "$NHP_MASK_NII" \
                    -expr "step(a - $THR)" \
                    -prefix "$TMPDIR/nhp_brain_mask_thr.nii.gz" \
                    -overwrite
                mv "$TMPDIR/nhp_brain_mask_thr.nii.gz" "$NHP_MASK_NII"

                # Erode slightly
                3dmask_tool \
                    -input "$NHP_MASK_NII" \
                    -prefix "$TMPDIR/nhp_brain_mask_eroded.nii.gz" \
                    -erode 1 \
                    -fill_holes \
                    -overwrite
                mv "$TMPDIR/nhp_brain_mask_eroded.nii.gz" "$NHP_MASK_NII"
                ;;

            deepbet)
                log " Running DeepBet..."
                deepbet \
                    -i "$TMPDIR/ecc_mean_3d.nii.gz" \
                    -o "$NHP_MASK_NII" \
                    -overwrite

                # DeepBet outputs a soft mask; threshold to binary
                log " Thresholding DeepBet output to binary mask..."
                local THR
                THR=$(3dmaskave -quiet -mask "$NHP_MASK_NII" "$NHP_MASK_NII" 2>/dev/null \
                      | awk '{print $1 * 0.5}')
                3dcalc \
                    -a "$NHP_MASK_NII" \
                    -expr "step(a - $THR)" \
                    -prefix "$TMPDIR/nhp_brain_mask_thr.nii.gz" \
                    -overwrite
                mv "$TMPDIR/nhp_brain_mask_thr.nii.gz" "$NHP_MASK_NII"

                # Erode slightly
                3dmask_tool \
                    -input "$NHP_MASK_NII" \
                    -prefix "$TMPDIR/nhp_brain_mask_eroded.nii.gz" \
                    -erode 1 \
                    -fill_holes \
                    -overwrite
                mv "$TMPDIR/nhp_brain_mask_eroded.nii.gz" "$NHP_MASK_NII"
                ;;
        esac

        # Verify mask dimensions match the ECC image
        MASK_VOX=$(mrinfo -size "$NHP_MASK_NII" | awk '{print $1, $2, $3}')
        ECC_VOX=$(mrinfo -size "$TMPDIR/ap_preproc.mif" | awk '{print $1, $2, $3}')
        log " NHP mask voxel dimensions (LxLyLz) : $MASK_VOX"
        log " ECC image voxel dimensions (LxLyLz) : $ECC_VOX"

        if [[ "$MASK_VOX" != "$ECC_VOX" ]]; then
            die "NHP mask dimensions ($MASK_VOX) do not match ECC image ($ECC_VOX). " \
                "The skull-stripping tool altered the geometry — check the method's output."
        fi

        # Convert to MIF for pipeline use and save to output
        mrconvert "$NHP_MASK_NII" "$TMPDIR/brain_mask.mif" -force
        mrconvert "$TMPDIR/brain_mask.mif" "$OUTDIR/${AP_BASE}_mask.nii.gz" -force

        log " Saved: ${AP_BASE}_mask.nii.gz (NHP auto-mask, method=$NHP_MASK_METHOD)"

    else

        # ------------------------------------------------------------------
        # DEFAULT: MRtrix3 dwi2mask (standard human / auto-masking path)
        # ------------------------------------------------------------------
        log "STEP 6: Creating brain mask from preprocessed DWI"

        dwi2mask "$TMPDIR/ap_preproc.mif" "$TMPDIR/brain_mask.mif" -force
        mrconvert "$TMPDIR/brain_mask.mif" "$OUTDIR/${AP_BASE}_mask.nii.gz" -force

        log " Saved: ${AP_BASE}_mask.nii.gz"

    fi

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
log "============================================================"