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


def clean_hdf_file(filename: Path):
    df = pandas.read_hdf(filename, "dynamics")
    df['count'] = df.groupby(['time', 'temperature', 'pressure'])['start_index'].transform('count')
    df = df.query('count > 10')
    df = df.query('time > 0')
    df = df.drop(columns=['count'], axis=1)

    df.to_hdf(filename.with_name(filename.stem + "_clean" + ".h5"), "dynamics")

    df = pandas.read_hdf(filename, "molecular_relaxations")
    df.to_hdf(filename.with_name(filename.stem + "_clean" + ".h5"), "molecular_relaxations")

if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = Path(sys.argv[1])
        clean_hdf_file(filename)
    else:
        print("Requires input path for processing.")


