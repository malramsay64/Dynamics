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
import sys
from pathlib import Path
from typing import Any, Tuple

import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
import click
import numpy as np
import pandas as pd
from sdanalysis.relaxation import series_relaxation_value


def _value(series: pd.Series):
    return bs.bootstrap(
        series.values, bs_stats.mean, alpha=0.5, num_iterations=1000
    ).value


def _lower(series: pd.Series):
    return bs.bootstrap(
        series.values, bs_stats.mean, alpha=0.5, num_iterations=1000
    ).lower_bound


def _upper(series: pd.Series):
    return bs.bootstrap(
        series.values, bs_stats.mean, alpha=0.1, num_iterations=1000
    ).upper_bound


def normalised_temperature(temperature: np.array, pressure: np.array) -> np.array:
    melting_points = {13.50: 1.35, "13.50": 1.35, 1.00: 0.36, "1.00": 0.36}
    temp_norm = np.full_like(temperature, np.nan)
    # Temperature can't be zero or below, so set these values to nan
    zero_mask = temperature > 0
    temp_norm[~zero_mask] = np.nan
    for p, t_m in melting_points.items():
        mask = np.logical_and(pressure == p, zero_mask)
        temp_norm[mask] = t_m / temperature[mask]

    return temp_norm


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
        count=df.groupby(["time", "temperature", "pressure"])["start_index"].transform(
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
    df_mol.index.names = ("keyframe", "molecule")
    df_mol = df_mol.reset_index()

    # Replace invalid values (2**32 - 1) with NaN's
    df_mol.replace(2 ** 32 - 1, np.nan, inplace=True)
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
        df.drop(columns="start_index")
        .groupby(["temperature", "pressure", "time"])
        .agg([_value, _lower, _upper])
    )
    df_agg.columns = ["".join(col).strip() for col in df_agg.columns.values]
    df_agg = df_agg.reset_index()
    df_agg["temp_norm"] = normalised_temperature(
        df_agg["temperature"], df_agg["pressure"]
    )

    df_agg.to_hdf(outfile, "dynamics")

    df_mol = pd.read_hdf(infile, "molecular_relaxations").drop(columns=["molecule"])

    # Taking the average over all the molecules from a single keyframe
    # makes the most sense from the perspective of computing an average,
    # since all molecules are present and not independent. One can be
    # fast because others are slow.
    df_mol = df_mol.groupby(["temperature", "pressure", "keyframe"]).mean()

    df_mol_agg = df_mol.groupby(["temperature", "pressure"]).agg(
        [_value, _lower, _upper]
    )
    df_mol_agg.columns = ["".join(col).strip() for col in df_mol_agg.columns.values]
    df_mol_agg = df_mol_agg.reset_index()
    df_mol_agg["temp_norm"] = normalised_temperature(
        df_mol_agg["temperature"], df_mol_agg["pressure"]
    )

    df_mol_agg.to_hdf(outfile, "molecular_relaxations")

    # Calculate the relaxation time from each keyframe
    df_relax = (
        df.set_index("time")
        .groupby(["temperature", "pressure", "start_index"])
        .agg(series_relaxation_value)
    )
    df_relax["inv_diffusion"] = 1 / df_relax["msd"]

    # Calculate the bootstrapped errors in the relaxation times
    df_relax_agg = df_relax.groupby(["temperature", "pressure"]).agg(
        [_value, _lower, _upper]
    )

    df_relax_agg.columns = ["".join(col).strip() for col in df_relax_agg.columns.values]
    df_relax_agg = df_relax_agg.reset_index()

    # Include temp_norm column.
    # This is the temperature normalised by the melting point
    df_relax_agg["temp_norm"] = normalised_temperature(
        df_relax_agg["temperature"], df_relax_agg["pressure"]
    )

    df_relax_agg.to_hdf(outfile, "relaxations")


@main.command()
@click.argument("output", type=click.Path(file_okay=True, dir_okay=False))
@click.argument(
    "infiles", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False)
)
def collate(output: Path, infiles: Tuple[Path, ...]) -> None:
    with pandas.HDFStore(output, "w") as dst:
        for file in infiles:
            with pandas.HDFStore(file) as src:
                for key in ["dynamics", "molecular_relaxations"]:
                    try:
                        dst.append(key, src.get(key))
                    except KeyError:
                        logger.warning("%s doesn't contain key %s", file, key)


if __name__ == "__main__":
    main()
