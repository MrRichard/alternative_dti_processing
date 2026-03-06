#!/usr/bin/env python3
"""
Create a native-space streamline density map from an ECC-corrected multi-shell DWI.

This script requires multi-shell data (>=2 non-zero b-value shells). If the input only
contains b=1000 (plus b0), it exits gracefully with a clear error message.
"""

import argparse
import os
import sys
from typing import List, Tuple

import nibabel as nib
import numpy as np
from dipy.core.gradients import gradient_table
from dipy.data import default_sphere
from dipy.direction import peaks_from_model
from dipy.io.gradients import read_bvals_bvecs
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel, auto_response_ssst
from dipy.segment.mask import median_otsu
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from dipy.tracking.streamline import Streamlines


def fail(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def info(msg: str) -> None:
    print(f"[INFO] {msg}")


def infer_gradient_files(ecc_path: str) -> Tuple[str, str]:
    base = ecc_path
    if base.endswith(".nii.gz"):
        base = base[:-7]
    elif base.endswith(".nii"):
        base = base[:-4]
    return f"{base}.bval", f"{base}.bvec"


def shell_centers(bvals: np.ndarray, tol: float = 100.0) -> List[float]:
    shells = sorted(float(b) for b in bvals if b > 50.0)
    if not shells:
        return []
    groups: List[List[float]] = [[shells[0]]]
    for b in shells[1:]:
        if abs(b - groups[-1][-1]) <= tol:
            groups[-1].append(b)
        else:
            groups.append([b])
    return [float(np.mean(group)) for group in groups]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate native-space streamline density map from ECC DWI."
    )
    parser.add_argument(
        "-i",
        "--ecc-dwi",
        required=True,
        help="ECC-corrected DWI NIfTI (for example, *_ECC.nii.gz).",
    )
    parser.add_argument(
        "-b",
        "--bval",
        default=None,
        help="bval file path. If omitted, inferred from ECC filename.",
    )
    parser.add_argument(
        "-v",
        "--bvec",
        default=None,
        help="bvec file path. If omitted, inferred from ECC filename.",
    )
    parser.add_argument(
        "-m",
        "--mask",
        default=None,
        help="Optional brain mask in DWI space. If omitted, estimated from b0.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output native-space streamline density map NIfTI.",
    )
    parser.add_argument(
        "--seed-density",
        type=int,
        default=1,
        help="Seeds per voxel edge for mask seeding grid (default: 1).",
    )
    parser.add_argument(
        "--step-size",
        type=float,
        default=0.5,
        help="Tracking step size in mm (default: 0.5).",
    )
    parser.add_argument(
        "--gfa-threshold",
        type=float,
        default=0.25,
        help="Stopping threshold on GFA map (default: 0.25).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not os.path.isfile(args.ecc_dwi):
        fail(f"ECC input not found: {args.ecc_dwi}")

    bval_path, bvec_path = args.bval, args.bvec
    if bval_path is None or bvec_path is None:
        inferred_bval, inferred_bvec = infer_gradient_files(args.ecc_dwi)
        bval_path = bval_path or inferred_bval
        bvec_path = bvec_path or inferred_bvec

    if not os.path.isfile(bval_path):
        fail(f"bval file not found: {bval_path}")
    if not os.path.isfile(bvec_path):
        fail(f"bvec file not found: {bvec_path}")

    info("Loading ECC DWI and gradient metadata")
    img = nib.load(args.ecc_dwi)
    data = img.get_fdata(dtype=np.float32)
    affine = img.affine

    if data.ndim != 4:
        fail("Input DWI must be 4D (x, y, z, volumes).")

    bvals, bvecs = read_bvals_bvecs(bval_path, bvec_path)
    if len(bvals) != data.shape[3]:
        fail(
            f"Volume count mismatch: DWI has {data.shape[3]} volumes, "
            f"but gradients define {len(bvals)} entries."
        )

    centers = shell_centers(bvals, tol=100.0)
    if len(centers) < 2:
        fail(
            "Input appears single-shell (one non-zero b-value shell). "
            "Multi-shell data is required for this post-processing step."
        )

    target_shell = max(centers)
    shell_window = 100.0
    keep_idx = np.where((bvals < 50.0) | (np.abs(bvals - target_shell) <= shell_window))[0]
    n_target = int(np.sum(np.abs(bvals - target_shell) <= shell_window))
    if n_target < 10:
        fail(
            f"Not enough diffusion directions in highest shell near b={target_shell:.1f}. "
            f"Found {n_target}, expected at least 10."
        )

    info(f"Detected non-zero shell centers: {', '.join(f'{c:.1f}' for c in centers)}")
    info(f"Using b0 + highest shell near b={target_shell:.1f} for CSD tracking")

    data_sub = data[..., keep_idx]
    bvals_sub = bvals[keep_idx]
    bvecs_sub = bvecs[keep_idx, :]
    gtab = gradient_table(bvals_sub, bvecs_sub, b0_threshold=50.0)

    if args.mask:
        if not os.path.isfile(args.mask):
            fail(f"Mask file not found: {args.mask}")
        mask_img = nib.load(args.mask)
        mask = mask_img.get_fdata() > 0
        if mask.shape != data.shape[:3]:
            fail(
                f"Mask shape {mask.shape} does not match DWI shape {data.shape[:3]}."
            )
        info("Using provided brain mask")
    else:
        b0_idx = np.where(bvals_sub < 50.0)[0]
        if b0_idx.size == 0:
            fail("No b=0 volumes available after shell selection.")
        _, mask = median_otsu(
            data_sub,
            vol_idx=b0_idx.tolist(),
            median_radius=3,
            numpass=2,
            autocrop=False,
            dilate=1,
        )
        info("Estimated brain mask using median_otsu on b0 volumes")

    info("Estimating response function and fitting CSD model")
    response, _ = auto_response_ssst(gtab, data_sub, roi_radii=10, fa_thr=0.7)
    if not np.isfinite(response[0]) or response[0] <= 0:
        fail("Could not estimate a valid CSD response function.")

    csd_model = ConstrainedSphericalDeconvModel(gtab, response)
    peaks = peaks_from_model(
        model=csd_model,
        data=data_sub,
        sphere=default_sphere,
        relative_peak_threshold=0.5,
        min_separation_angle=25,
        mask=mask,
        return_sh=True,
        parallel=False,
    )

    stopping = ThresholdStoppingCriterion(peaks.gfa, args.gfa_threshold)
    seeds = utils.seeds_from_mask(mask, affine=affine, density=args.seed_density)
    if len(seeds) == 0:
        fail("No seeds were generated from the mask.")

    info("Running deterministic local tracking")
    streamline_generator = LocalTracking(
        direction_getter=peaks,
        stopping_criterion=stopping,
        seeds=seeds,
        affine=affine,
        step_size=args.step_size,
    )
    streamlines = Streamlines(streamline_generator)
    if len(streamlines) == 0:
        fail("Tracking produced zero streamlines.")

    info(f"Generated {len(streamlines)} streamlines")
    density = utils.density_map(streamlines, affine=affine, vol_dims=data.shape[:3])
    density_img = nib.Nifti1Image(density.astype(np.float32), affine, img.header)
    nib.save(density_img, args.output)
    info(f"Wrote native-space streamline density map: {args.output}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        fail("Interrupted by user.", code=130)
