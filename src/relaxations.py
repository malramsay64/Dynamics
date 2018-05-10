#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.
"""Compute a single relaxation value from the dynamics."""

import numpy as np
import pandas as pd
from scipy.stats import hmean
from sdanalysis import relaxation


def max_relaxation(time: np.ndarray, value: np.ndarray) -> float:
    """Time at which the maximum value is recorded.

    Args:
        time (np.ndarray): The time index
        value (np.ndarray): The value at each of the time indices

    Returns:
        float: The time at which the maximum value occurs.
        float: Value of the maximum.

    """
    max_val_index = np.argmax(value)
    return time[max_val_index], value[max_val_index]


def compute_relaxations(series: pd.Series):
    index = series.index.get_level_values(0)
    if series.name in ['msd']:
        relax, relax_error = relaxation.diffusion_constant(index, series.values)
    elif series.name in ['struct_msd']:
        relax, relax_error = relaxation.threshold_relaxation(
            index, series.values, threshold=0.16, greater=False
        )
    elif series.name in ['alpha', 'gamma']:
        relax, relax_error = max_relaxation(index, series.values)
    else:
        relax, err_max, err_min = relaxation.exponential_relaxation(
            index, series.values
        )
        relax_error = err_max - err_min
    return relax


def translate_relaxation(quantity: str) -> str:
    translation = {
        'alpha': 'max_alpha_time',
        'gamma': 'max_gamma_time',
        'com_struct': 'tau_F',
        'msd': 'diffusion_constant',
        'rot1': 'tau_R1',
        'rot2': 'tau_R2',
        'struct': 'tau_S',
    }
    return translation.get(quantity, quantity)


def main():
    df_dyn = pd.read_hdf('data/analysis/dynamics.h5', 'dynamics')
    # Remove columns with no relaxation value to calculate
    df_dyn.drop(
        ['mean_displacement', 'mean_rotation', 'mfd', 'overlap', 'start_index'],
        inplace=True,
    )
    # Average over all intial times
    df_dyn = df_dyn.groupby(['time', 'temperature', 'pressure']).mean()
    relaxations = df_dyn.groupby(['temperature', 'pressure']).aggregate(
        compute_relaxations
    )
    relaxations.columns = [
        translate_relaxation(quantity) for quantity in relaxations.columns
    ]
    df_mol = pd.read_hdf('data/analysis/dynamics.h5', 'molecular_relaxations')
    df_mol.replace(2 ** 32 - 1, np.nan, inplace=True)
    df_mol.index.names = ['init_frame', 'molecule']
    df_mol = df_mol.groupby(['init_frame', 'temperature', 'pressure']).mean(
        skipna=False
    )
    df_mol = df_mol.groupby(['temperature', 'pressure']).agg(['mean', hmean])
    df_mol.columns = ['_'.join(f) for f in df_mol.columns.tolist()]
    pd.concat([df_mol, relaxations], axis=1).to_hdf(
        'data/analysis/dynamics.h5', 'relaxations'
    )


if __name__ == '__main__':
    main()
