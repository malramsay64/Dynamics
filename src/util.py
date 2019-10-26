#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Utility functions for analysis across the rest of the files."""

from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


def _read_temperatures(filename: Path) -> Dict[float, float]:
    """Read temperatures from a CSV file and format for simple translation.

    Args:
        filename: An input file which contains the melting points for each pressure.

    """
    df = pd.read_csv(filename)
    melting_points = {}
    for _, row in df.iterrows():
        melting_points[float(row["pressure"])] = float(row["melting_point"])
    return melting_points


def normalised_temperature(temperature: np.array, pressure: np.array) -> np.array:
    """Calculate a normalised temperature based on the melting point.

    This calculates normalised temperatures based on the melting points specified in the
    'results/melting_points.csv' file.

    Args:
        temperature: A collection of temperatures for which the normalised 
            value is to be calculated. These are all expected to be floating
            point values.
        pressure: The pressure corresponding to each of the input temperatures.
            These are all expected to be floating point values.

    Returns:
        Array containing the temperatures normalised by the appropriate melting
        point for the given pressure.

    """
    # Ensure input is floating point values.
    if not np.issubdtype(temperature.dtype, np.floating):
        raise ValueError(
            "The temperature needs to be specified as floating point values."
        )
    if not np.issubdtype(pressure.dtype, np.floating):
        raise ValueError("The pressure needs to be specified as floating point values.")

    melting_points = _read_temperatures(
        Path(__file__).parent.parent / "results/melting_points.csv"
    )

    # Initialise output array
    temp_norm = np.full_like(temperature, np.nan)

    # Temperature can't be zero or below, so set these values to nan
    zero_mask = temperature > 0
    temp_norm[~zero_mask] = np.nan

    # Iterate over all pressures in melting_points dictionary
    for p, t_m in melting_points.items():
        # Find values with the same pressure
        mask = np.logical_and(pressure == p, zero_mask)
        # Calculate normalised temperature
        temp_norm[mask] = temperature[mask] / t_m

    return temp_norm
