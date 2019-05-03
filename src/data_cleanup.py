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
import pandas
import sys
from pathlib import Path
import click


@click.command()
@click.argument("infile", type=click.Path(file_okay=True, dir_okay=False, exists=True))
@click.option(
    "--min-samples",
    default=10,
    type=int,
    help="Minimum number of samples for each data point.",
)
def clean_hdf_file(infile: Path, min_samples: int):
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

    df.to_hdf(infile.with_name(infile.stem + "_clean" + ".h5"), "dynamics")

    df = pandas.read_hdf(infile, "molecular_relaxations")
    df.to_hdf(infile.with_name(infile.stem + "_clean" + ".h5"), "molecular_relaxations")


if __name__ == "__main__":
    clean_hdf_file()
