#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from src import calc_dynamics


@pytest.mark.parametrize("press, temp", [("13.50", 1.35), ("1.00", 0.36) , (13.50, 1.35), (1.00, 0.36)])
def test_temperature_normalisation_melting(press, temp):
    shape = 10
    pressure = np.full(shape, press)
    temperature = np.full(shape, temp)
    result = calc_dynamics.normalised_temperature(temperature, pressure)
    assert_array_almost_equal(result, np.ones(shape))


@pytest.mark.parametrize("press, temp", [("13.50", 1.35), ("1.00", 0.36) , (13.50, 1.35), (1.00, 0.36)])
def test_temperature_normalisation_one(press, temp):
    shape = 10
    pressure = np.full(shape, press)
    temperature = np.ones(shape)
    result = calc_dynamics.normalised_temperature(temperature, pressure)
    assert_array_almost_equal(result, np.full(shape, temp))


@pytest.mark.parametrize("press, temp", [("13.50", 1.35), ("1.00", 0.36) , (13.50, 1.35), (1.00, 0.36)])
def test_temperature_normalisation_zero(press, temp):
    shape = 10
    pressure = np.full(shape, press)
    temperature = np.zeros(shape)
    result = calc_dynamics.normalised_temperature(temperature, pressure)
    assert_array_almost_equal(result, np.full(shape, np.nan))
