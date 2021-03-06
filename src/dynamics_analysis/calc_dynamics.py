#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Clean the raw quantities from the dynamics analysis.

This is a collection of utilities for cleaning up the raw data from the calculation of
the dynamics.

"""
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import click
import numpy as np
import pandas as pd
import scipy.stats
import sdanalysis
from sdanalysis.relaxation import series_relaxation_value

from .util import normalised_temperature

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def _read_temperatures(filename: Path) -> Dict[float, float]:
    """Read temperatures from a CSV file and format for simple translation.

    Args:
        filename: An input file which contains the melting points for each pressure.

    """
    df = pd.read_csv(filename)
    melting_points = {}
    for row in df:
        melting_points[float(row["pressure"])] = float(row["melting_point"])
    return melting_points


def inv_mean(values: pd.Series) -> float:
    return (1 / values).mean()


def inv_sem(values: pd.Series) -> float:
    return (1 / values).sem()


@click.group()
def main():
    pass


@main.command()
@click.argument("infile", type=click.Path(file_okay=True, dir_okay=False, exists=True))
@click.option(
    "--min-samples",
    default=10,
    type=int,
    help="Minimum number of samples for each data point.",
)
def clean(infile: Path, min_samples: int):
    infile = Path(infile)
    # Cleanup the dynamics dataset
    df = pd.read_hdf(infile, "dynamics")

    # Most of the values are plotted on a log scale for the time axis, values less than
    # or equal to 0 cause issues.
    df = df.query("time > 0")

    # We want to discard values where there are not enough to get decent statistics, in
    # this case I have chosen 10 as the magic number.
    df = df.assign(
        count=df.groupby(["time", "temperature", "pressure"])["keyframe"].transform(
            "count"
        )
    )
    df = df.query("count > @min_samples")
    #  # Don't want the count in the final dataset, just a temporary column
    df = df.drop(columns=["count"], axis=1)

    # The values where the MSD is greater than 100 are going to have issues with
    # the periodic boundary conditions so remove those columns.
    df = df.query("msd < 100")
    df = df.reset_index()

    df.to_hdf(infile.with_name(infile.stem + "_clean" + ".h5"), "dynamics")

    df_mol = pd.read_hdf(infile, "molecular_relaxations")
    df_mol = df_mol.reset_index()

    # Replace invalid values (2**32 - 1) with NaN's
    df_mol.replace((2 ** 32 - 1) * 0.005, np.nan, inplace=True)
    # Remove keyframes where relaxation hasn't completed,
    # that is there are NaN values present.
    df_mol = df_mol.groupby(["keyframe", "temperature", "pressure"]).filter(
        lambda x: x.isna().sum().sum() == 0
    )

    df_mol = df_mol.assign(
        count=df_mol.groupby(["temperature", "pressure"])["keyframe"].transform("count")
    )
    df_mol = df_mol.query("count > @min_samples")
    #  # Don't want the count in the final dataset, just a temporary column
    df_mol = df_mol.drop(columns=["count"], axis=1)

    df_mol.to_hdf(
        infile.with_name(infile.stem + "_clean" + ".h5"), "molecular_relaxations"
    )


@main.command()
@click.argument("infile", type=click.Path(file_okay=True, dir_okay=False, exists=True))
def bootstrap(infile):
    infile = Path(infile)
    outfile = infile.with_name(infile.stem + "_agg" + ".h5")
    df = pd.read_hdf(infile, "dynamics").drop(columns=["index"])

    df_agg = (
        df.drop(columns="keyframe")
        .groupby(["temperature", "pressure", "time"])
        .agg(["mean", "sem"])
    )
    df_agg.columns = ["_".join(col).strip() for col in df_agg.columns.values]
    df_agg = df_agg.reset_index()
    df_agg["inv_temp_norm"] = 1 / normalised_temperature(
        df_agg["temperature"].values, df_agg["pressure"].values
    )

    df_agg.to_hdf(outfile, "dynamics")

    df_mol = pd.read_hdf(infile, "molecular_relaxations").drop(columns=["molecule"])

    df_mol = df_mol.drop(columns=["keyframe", "index"])
    df_mol_agg = df_mol.groupby(["temperature", "pressure"]).agg(
        ["mean", "sem", inv_mean, inv_sem]
    )

    df_mol_agg.columns = ["_".join(col).strip() for col in df_mol_agg.columns.values]
    df_mol_agg = df_mol_agg.reset_index()
    df_mol_agg["inv_temp_norm"] = 1 / normalised_temperature(
        df_mol_agg["temperature"].values, df_mol_agg["pressure"].values
    )

    df_mol_agg.to_hdf(outfile, "molecular_relaxations")

    df["msr"] = df["msr"].mask(lambda x: x < 10, 1000)
    # Calculate the relaxation time from each keyframe
    df_relax = (
        df.set_index("time")
        .groupby(["temperature", "pressure", "keyframe"])
        .agg(series_relaxation_value)
    )
    df_relax["inv_diffusion"] = 1 / df_relax["msd"]
    df_relax["inv_diffusion_rot"] = 1 / df_relax["msr"]

    # Calculate the bootstrapped errors in the relaxation times
    df_relax_agg = df_relax.groupby(["temperature", "pressure"]).agg(["mean", "sem"])

    df_relax_agg.columns = [
        "_".join(col).strip() for col in df_relax_agg.columns.values
    ]
    df_relax_agg = df_relax_agg.reset_index()

    # Include inv_temp_norm column.
    # This is the temperature normalised by the melting point
    df_relax_agg["inv_temp_norm"] = 1.0 / normalised_temperature(
        df_relax_agg["temperature"].values, df_relax_agg["pressure"].values
    )

    df_relax_agg.to_hdf(outfile, "relaxations")


@main.command()
@click.argument("output", type=click.Path(file_okay=True, dir_okay=False))
@click.argument(
    "infiles", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False)
)
def collate(output: Path, infiles: Tuple[Path, ...]) -> None:
    with pd.HDFStore(output, "w") as dst:
        for file in infiles:
            with pd.HDFStore(file) as src:
                for key in ["dynamics", "molecular_relaxations"]:
                    try:
                        df = src.get(key)
                    except KeyError:
                        logger.warning("File %s doesn't contain key %s", file, key)

                    if key == "dynamics":
                        # The timestep is given the time column, so convert that here
                        df["timestep"] = df["time"]
                        # This converts the timestep to the real time
                        df["time"] = df["timestep"] * 0.005
                    elif key == "molecular_relaxations":
                        df = df.set_index(
                            ["pressure", "temperature", "keyframe", "molecule"]
                        )
                        df *= 0.005
                        df = df.reset_index()

                    df["temperature"] = df["temperature"].astype(float)
                    df["pressure"] = df["pressure"].astype(float)
                    dst.append(key, df)


@main.command()
@click.argument("infile", type=click.Path(file_okay=True, exists=True, dir_okay=False))
@click.argument("outfile", type=click.Path(file_okay=True, dir_okay=False))
def stokes_einstein(infile: Path, outfile: Path):
    sdanalysis.read.process_file(infile, outfile=outfile, wave_number=2.90)

    data: List[Dict[str, np.ndarray]] = []
    dyn = None
    for frame in sdanalysis.read.open_trajectory(infile, progressbar=True):
        if dyn is None:
            dyn = sdanalysis.dynamics.Dynamics.from_frame(frame)
        data.append(
            pd.DataFrame(
                {
                    "time": dyn.compute_time_delta(frame.timestep) * 0.005,
                    "timestep": dyn.compute_time_delta(frame.timestep),
                    "molecule": dyn.get_molid(),
                    "displacement": dyn.get_displacements(frame.position),
                    "rotation": dyn.get_rotations(frame.orientation),
                }
            )
        )
    disp_df = pd.concat(data)
    disp_df.to_hdf(outfile, "displacements")


if __name__ == "__main__":
    main()
