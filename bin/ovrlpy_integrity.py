#!/usr/bin/env python3
"""Run ovrlpy signal-integrity QC on a Xenium bundle.

Vertical signal integrity (VSI) is low where cells overlap in z, so no 2D
segmentation can separate their transcripts. Writes into --outdir:

    signal_integrity.png   VSI map with histogram
    integrity_map.npy      VSI per grid pixel
    signal_map.npy         transcript signal per grid pixel
    grid.json              origin and pixel size (um) of both maps
    doublets.csv           local VSI minima: likely overlapping cells
    cell_integrity.csv     per-cell VSI and signal, using the bundle's cell assignment
"""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import ovrlpy
import polars as pl

matplotlib.use("Agg")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bundle", type=Path, required=True, help="Xenium output bundle")
    parser.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
    parser.add_argument("--min-qv", type=float, default=20, help="Minimum transcript Q-score")
    parser.add_argument("--threads", type=int, default=1, help="Worker threads")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    return parser.parse_args()


def cell_integrity(ovrlp: ovrlpy.Ovrlp) -> pl.DataFrame:
    """Summarise VSI and signal per cell from the transcripts' cell assignment.

    Args:
        ovrlp: Fitted ovrlpy model whose transcripts keep the ``cell_id`` column.

    Returns:
        One row per cell: mean and minimum VSI, mean signal, pixels covered.
    """
    pixels = ovrlpy.cell_integrity_from_transcripts(ovrlp, cell_id="cell_id", unassigned="UNASSIGNED")
    return (
        pixels.group_by("cell_id")
        .agg(
            vsi_mean=pl.col("vsi").mean(),
            vsi_min=pl.col("vsi").min(),
            signal_mean=pl.col("signal").mean(),
            n_pixels=pl.len(),
        )
        .sort("cell_id")
    )


def main() -> None:
    """Fit ovrlpy and write maps, plot, doublet calls and per-cell integrity."""
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    transcripts = ovrlpy.io.read_Xenium(
        args.bundle / "transcripts.parquet",
        min_qv=args.min_qv,
        additional_columns=["cell_id"],
        n_threads=args.threads,
    ).with_columns(pl.col("cell_id").cast(pl.String))

    ovrlp = ovrlpy.Ovrlp(transcripts, n_workers=args.threads, random_state=args.seed)
    ovrlp.analyse()

    ovrlpy.plot_signal_integrity(ovrlp).savefig(args.outdir / "signal_integrity.png", dpi=200)
    np.save(args.outdir / "integrity_map.npy", ovrlp.integrity_map)
    np.save(args.outdir / "signal_map.npy", ovrlp.signal_map)
    grid = {"origin": [float(v) for v in ovrlp.origin], "gridsize": float(ovrlp.gridsize)}
    (args.outdir / "grid.json").write_text(json.dumps(grid, indent=2))
    ovrlp.detect_doublets().write_csv(args.outdir / "doublets.csv")
    cell_integrity(ovrlp).write_csv(args.outdir / "cell_integrity.csv")


if __name__ == "__main__":
    main()
