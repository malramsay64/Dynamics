#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

import numpy as np

from ..relaxations import max_relaxation


def test_max_relaxation():
    time = np.arange(101)
    value = 50 - np.abs(np.arange(-50, 50))
    max_time, error = max_relaxation(time, value)
    assert max_time == 50
    assert error == 1

def test_max_relaxation_nan():
    time = np.arange(101)
    value = 50. - np.abs(np.arange(-50, 50))
    value[0] = np.nan
    max_time, error = max_relaxation(time, value)
    assert max_time == 50
    assert error == 1
