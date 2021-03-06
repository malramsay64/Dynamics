#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""An implementation of brownian motion using stochastic processes.

This is a simulation of Brownian motion, allowing for the calculation of quantities for
pure brownian motion for comparison to the dynamics quantities I am calculating for my
trimer simulations.

The Brownian motion simulation is taken from the scipy-cookbook (https://scipy-cookbook.readthedocs.io/items/BrownianMotion.html).

"""

from math import sqrt

import numpy as np
from scipy.stats import norm
from sdanalysis import dynamics


def get_brownian_relax(
    steps: int = 1000, time: int = 10, step_size: float = 0.25, molecules: int = 2000
):
    """Compute the translational relaxation times for Brownian motion."""
    # Change in position from motion
    delta = np.linalg.norm(
        brownian(np.zeros((2, molecules)), steps, time / steps, step_size), axis=0,
    )
    # These are the relaxation quantities to calculate
    tau_F = dynamics.MolecularRelaxation(molecules, 0.4)
    tau_D = dynamics.MolecularRelaxation(molecules, 2.0)
    tau_L = dynamics.LastMolecularRelaxation(molecules, 0.4, 2.0)
    # Add displacements to relaxation calculations
    for time in range(delta.shape[-1]):
        tau_F.add(time, delta[:, time])
        tau_D.add(time, delta[:, time])
        tau_L.add(time, delta[:, time])
    return tau_F, tau_D, tau_L


def brownian(x0: np.ndarray, n: int, dt: float, delta: float, out=None):
    """Generate an instance of Brownian motion (i.e. the Wiener process):

        X(t) = X(0) + N(0, delta**2 * t; 0, t)

    where N(a,b; t0, t1) is a normally distributed random variable with mean a and
    variance b.  The parameters t0 and t1 make explicit the statistical
    independence of N on different time intervals; that is, if [t0, t1) and
    [t2, t3) are disjoint intervals, then N(a, b; t0, t1) and N(a, b; t2, t3)
    are independent.

    Written as an iteration scheme,

        X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)


    If `x0` is an array (or array-like), each value in `x0` is treated as
    an initial condition, and the value returned is a numpy array with one
    more dimension than `x0`.

    Arguments
    ---------
    x0 : float or numpy array (or something that can be converted to a numpy array
         using numpy.asarray(x0)).
        The initial condition(s) (i.e. position(s)) of the Brownian motion.
    n : int
        The number of steps to take.
    dt : float
        The time step.
    delta : float
        delta determines the "speed" of the Brownian motion.  The random variable
        of the position at time t, X(t), has a normal distribution whose mean is
        the position at time t=0 and whose variance is delta**2*t.
    out : numpy array or None
        If `out` is not None, it specifies the array in which to put the
        result.  If `out` is None, a new numpy array is created and returned.

    Returns
    -------
    A numpy array of floats with shape `x0.shape + (n,)`.

    Note that the initial value `x0` is not included in the returned array.
    """

    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    r = norm.rvs(size=x0.shape + (n,), scale=delta * sqrt(dt))

    # If `out` was not given, create an output array.
    if out is None:
        out = np.empty(r.shape)

    # This computes the Brownian motion by forming the cumulative sum of
    # the random samples.
    np.cumsum(r, axis=-1, out=out)

    # Add the initial condition.
    out += np.expand_dims(x0, axis=-1)

    return out


def brownian_rotation(x0: np.ndarray, n: int, dt: float, delta: float, out=None):
    """Generate an instance of rotational Brownian motion (i.e. the Wiener process):

        X(t) = X(0) + N(0, delta**2 * t; 0, t)

    where N(a,b; t0, t1) is a normally distributed random variable with mean a and
    variance b.  The parameters t0 and t1 make explicit the statistical
    independence of N on different time intervals; that is, if [t0, t1) and
    [t2, t3) are disjoint intervals, then N(a, b; t0, t1) and N(a, b; t2, t3)
    are independent.

    Written as an iteration scheme,

        X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)


    If `x0` is an array (or array-like), each value in `x0` is treated as
    an initial condition, and the value returned is a numpy array with one
    more dimension than `x0`.

    This is different from the brownian function in that the values are 
    restricted to the range :math:`[0, 2\pi)`.

    Arguments
    ---------
    x0 : float or numpy array (or something that can be converted to a numpy array
         using numpy.asarray(x0)).
        The initial condition(s) (i.e. position(s)) of the Brownian motion.
    n : int
        The number of steps to take.
    dt : float
        The time step.
    delta : float
        delta determines the "speed" of the Brownian motion.  The random variable
        of the position at time t, X(t), has a normal distribution whose mean is
        the position at time t=0 and whose variance is delta**2*t.
    out : numpy array or None
        If `out` is not None, it specifies the array in which to put the
        result.  If `out` is None, a new numpy array is created and returned.

    Returns
    -------
    A numpy array of floats with shape `x0.shape + (n,)`.

    Note that the initial value `x0` is not included in the returned array.
    """

    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    r = norm.rvs(size=x0.shape + (n,), scale=delta * sqrt(dt))

    # Convert each element to a complex number
    r = np.exp(1j * r)

    # If `out` was not given, create an output array.
    if out is None:
        out = np.empty_like(r)

    # This computes the Brownian motion by forming the cumulative product of
    # the random angular jumps.
    np.cumprod(r, axis=-1, out=out)

    # Add the initial condition.
    out *= np.expand_dims(x0, axis=-1)

    return out
