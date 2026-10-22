#!/usr/bin/env python3
"""Segment a Xenium morphology image with InstanSeg.

Writes nuclei and cell label masks (uint32 TIFF, full morphology resolution,
pixel units): the format Xenium Ranger import-segmentation expects, matching
how nf-core/spatialaxe imports Cellpose masks.

Inference runs on a pyramid level (default 1, i.e. half resolution, close to
the 0.5 um pixel size InstanSeg models are trained at) to bound memory; masks
are upscaled back to full resolution with nearest-neighbour repetition.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import tifffile
import torch
from instanseg import InstanSeg


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bundle", type=Path, required=True, help="Xenium output bundle")
    parser.add_argument("--image", type=Path, help="Morphology image to use instead of the bundle's")
    parser.add_argument("--model", default="fluorescence_nuclei_and_cells", help="InstanSeg model")
    parser.add_argument("--level", type=int, default=1, help="Pyramid level for inference (0 = full)")
    parser.add_argument("--channels", type=int, nargs="+", help="Channel indices to use (default: all)")
    parser.add_argument("--nuclei", type=Path, default=Path("nuclei.tif"), help="Output nuclei mask")
    parser.add_argument("--cells", type=Path, default=Path("cells.tif"), help="Output cell mask")
    return parser.parse_args()


def morphology_files(bundle: Path, image: Path | None) -> list[Path]:
    """Find the 2D morphology image file(s), one channel per file.

    Uses the same precedence as nf-core/spatialaxe: an explicit image, then
    ``morphology_focus/`` (XOA >= 2.0, one file per stain), then
    ``morphology_focus.ome.tif`` (XOA 1.x, DAPI only).

    Args:
        bundle: Xenium output bundle.
        image: Optional explicit image path.

    Returns:
        Image files, in channel order.
    """
    if image is not None:
        return [image]
    focus_dir = bundle / "morphology_focus"
    if focus_dir.is_dir():
        return sorted(focus_dir.glob("*.ome.tif"))
    if (bundle / "morphology_focus.ome.tif").exists():
        return [bundle / "morphology_focus.ome.tif"]
    raise SystemExit(f"No 2D morphology_focus image in {bundle}; pass --image")


def read_channels(files: list[Path], level: int) -> tuple[np.ndarray, tuple[int, int]]:
    """Read every channel at one pyramid level.

    Files are read without OME multi-file resolution so each file yields only
    its own channel(s).

    Args:
        files: Image files from :func:`morphology_files`.
        level: Pyramid level to read.

    Returns:
        Image stack (channels, y, x) as float32, and the full-resolution (y, x) shape.
    """
    planes, full_shape = [], None
    for f in files:
        with tifffile.TiffFile(f, is_ome=False) as tif:
            series = tif.series[0]
            full_shape = full_shape or tuple(series.levels[0].shape[-2:])
            img = series.levels[min(level, len(series.levels) - 1)].asarray()
        planes.append(img.reshape(-1, *img.shape[-2:]))
    return np.concatenate(planes).astype(np.float32), full_shape


def upscale(mask: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Nearest-neighbour upscale a label mask to ``shape``.

    Args:
        mask: Label mask (y, x) at a lower pyramid level.
        shape: Full-resolution (y, x) shape.

    Returns:
        uint32 label mask of exactly ``shape``.
    """
    factor = int(round(shape[0] / mask.shape[0]))
    if factor > 1:
        mask = np.repeat(np.repeat(mask, factor, axis=0), factor, axis=1)
    out = np.zeros(shape, dtype=np.uint32)
    h, w = min(shape[0], mask.shape[0]), min(shape[1], mask.shape[1])
    out[:h, :w] = mask[:h, :w]
    return out


def main() -> None:
    """Run InstanSeg and write nuclei and cell masks."""
    args = parse_args()
    pixel_size = json.loads((args.bundle / "experiment.xenium").read_text())["pixel_size"]
    files = morphology_files(args.bundle, args.image)
    image, full_shape = read_channels(files, args.level)
    if args.channels:
        image = image[args.channels]
    level_pixel_size = pixel_size * full_shape[0] / image.shape[1]
    print(f"{len(files)} file(s), image {image.shape} at {level_pixel_size:.4f} um/px")

    model = InstanSeg(args.model, verbosity=1)
    with torch.no_grad():
        labels = model.eval_medium_image(
            image, pixel_size=level_pixel_size, target="all_outputs", return_image_tensor=False
        )
    labels = labels.cpu().numpy().reshape(-1, *image.shape[1:])  # (nuclei, cells, y, x)

    for name, mask, path in zip(("nuclei", "cells"), labels, (args.nuclei, args.cells)):
        mask = upscale(mask, full_shape)
        tifffile.imwrite(path, mask, compression="zlib", bigtiff=True)
        print(f"{name}: {len(np.unique(mask)) - 1} objects -> {path}")


if __name__ == "__main__":
    main()
