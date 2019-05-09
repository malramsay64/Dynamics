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
from typing import Any

import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
import click
import numpy as np
import pandas

from sdanalysis.relaxation import series_relaxation_value


def _value(series: pandas.Series):
    return bs.bootstrap(
        series.values, bs_stats.mean, alpha=0.5, num_iterations=1000
    ).value


def _lower(series: pandas.Series):
    return bs.bootstrap(
        series.values, bs_stats.mean, alpha=0.5, num_iterations=1000
    ).lower_bound


def _upper(series: pandas.Series):
    return bs.bootstrap(
        series.values, bs_stats.mean, alpha=0.1, num_iterations=1000
    ).upper_bound


def bootstap_series(
    series: pandas.Series,
    stat_function: Any,
    alpha: float = 0.5,
    num_iterations: int = 1000,
) -> bs.BootstrapResults:
    results = np.empty(num_iterations)
    num_samples = len(series)
    groups = series.groupby(level=0)
    keys = list(groups.groups.keys())
    for i in range(num_iterations):
        results[i] = stat_function(
            pandas.concat(
                groups.get_group(keys[key])
                for key in np.random.choice(len(keys), len(keys), replace=True)
            )
        )

    low = np.percentile(results, 100 * (alpha / 2.0))
    val = np.percentile(results, 50)
    high = np.percentile(results, 100 * (1 - alpha / 2.0))

    return bs.BootstrapResults(low, val, high)


def _relax_value(series: pandas.Series):
    return bootstap_series(
        series, series_relaxation_value, alpha=0.1, num_iterations=1000
    ).value


def _relax_lower(series: pandas.Series):
    return bootstap_series(
        series, series_relaxation_value, alpha=0.1, num_iterations=1000
    ).lower_bound


def _relax_upper(series: pandas.Series):
    return bootstap_series(
        series, series_relaxation_value, alpha=0.1, num_iterations=1000
    ).upper_bound


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
    df = pandas.read_hdf(infile, "dynamics")

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

    df_mol = pandas.read_hdf(infile, "molecular_relaxations")
    df_mol.index.names = ("keyframe", "molecule")
    df_mol = df_mol.reset_index()

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
    df = pandas.read_hdf(infile, "dynamics")

    df_agg = df.groupby(["temperature", "pressure", "time"]).agg(
        [_value, _lower, _upper]
    )
    df_agg.columns = ["".join(col).strip() for col in df_agg.columns.values]
    df_agg = df_agg.reset_index()

    df_agg.to_hdf(outfile, "dynamics")

    df_mol = pandas.read_hdf(infile, "molecular_relaxations")
    df_mol = df_mol.reset_index()

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

    df_mol_agg.to_hdf(outfile, "molecular_relaxations")

    df_relax = df.groupby(["temperature", "pressure", "start_index"]).agg(
        series_relaxation_value
    )
    df_relax["inv_diffusion"] = 1 / df_relax["msd"]
    df_relax_agg = df_relax.groupby(["temperature", "pressure"]).agg(
        [_value, _lower, _upper]
    )

    df_relax_agg.columns = ["".join(col).strip() for col in df_relax_agg.columns.values]
    df_relax_agg = df_relax_agg.reset_index()

    df_relax_agg.to_hdf(outfile, "relaxations")


if __name__ == "__main__":
    main()
