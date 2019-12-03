#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Model the Vogel-Tamman-Fulcher relations for relaxation timescales.

This module provides functions for modelling the Vogel-Tamman-Fulcher equation and fitting results from
experimental systems, enabling the extrapolation of the glass transition temperature and the fragility.
"""

from typing import Optional, Tuple

import numpy as np
import scipy.misc
import scipy.optimize


def vogel_tamman_fulcher(
    inv_T_norm: np.ndarray, tau0: float, D: float, T0: float
) -> np.ndarray:
    r"""Function describing the Vogel-Tamman-Fulcher relations.

    This is a relationship which describes the dynamics as they approach the 
    glass transition temperature. Once fit it allows the prediction of
    the relaxation time :math:`\tau_\alpha`.

    .. math::

        \tau_\alpha = \tau_0 \exp\left[\frac{DT_0}{T-T_0}\right] 

    Args:
        inv_T_norm: The normalised temperature for which to find the values.
        tau0: This is an extrapolated relaxation time found through fitting to
            calculated data.
        D: A breadth parameter found through fitting to a dataset.
        T0: An extrapolated temperature found through fitting.

    Returns:
        The predicted values based on the input.

    """
    return tau0 * np.exp(D * T0 / (1 / inv_T_norm - T0))


def fit_vtf(
    inv_temp_norm: np.ndarray, values: np.ndarray, errors: Optional[np.ndarray] = None
):
    """Fit the Vogel-Tamman-Fulcher relations for a specific relaxation quantity.

    Args:
        inv_temp_norm: The normalised temperature.
        quantity: The quantity to find the relation. This will be without the '_mean' suffix.

    """
    fit, error = scipy.optimize.curve_fit(
        vogel_tamman_fulcher, inv_temp_norm, values, sigma=errors, p0=(1.0, 1.0, 0.0),
    )
    return fit, np.sqrt(np.diag(error))


def find_glass_transition_fragility(
    inv_temp_norm: np.ndarray,
    values: np.ndarray,
    errors: Optional[np.ndarray] = None,
    root=1e14,
    bracket=(1.2, 1.4),
) -> Tuple[float, float]:
    fit, _ = fit_vtf(inv_temp_norm, values, errors)
    root = scipy.optimize.root_scalar(
        lambda x: vogel_tamman_fulcher(x, *fit) - root, bracket=bracket,
    ).root
    fragility = scipy.misc.derivative(
        lambda x: np.log(vogel_tamman_fulcher(x, *fit)), root, dx=1e-6
    )
    return root, fragility
