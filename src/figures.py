#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Code to assist in the creation of figures."""

import logging
from typing import Optional

import altair as alt
import click
import freud
import matplotlib.pyplot as plt
import numpy as np

import sdanalysis

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def my_theme() -> Dict[str, Any]:
    """Define an Altair theme to use in all visualisations.

    This defines a simple theme, specifying the aspect ratio of 4:6 
    and removing the grid from the figure which is distracting.

    """
    return {"config": {"view": {"height": 400, "width": 600}, "axis": {"grid": False}}}


def use_my_theme(alt):
    """Register and my custom Altair theme."""
    # register and enable the theme
    alt.themes.register("my_theme", my_theme())
    alt.themes.enable("my_theme")


def radial_distribution(
    infile,
    dr: Optional[float] = None,
    rmax: Optional[float] = None,
    num_frames: int = 1000,
) -> freud.density.RDF:
    if rmax is None:
        snapshot = next(sdanalysis.open_trajectory(infile))
        box = snapshot.box
        # 2D Box
        if box[2] == 1.0:
            rmax = np.min(box[:2]) / 2.2
        # 3D box
        else:
            rmax = np.min(box[:3]) / 2.2
        logger.debug("rmax set to None, calculated best value as %s", rmax)

    if dr is None:
        dr = rmax / 1000
        logger.debug("dr set to None, calculated best value as %s", dr)

    logger.debug("rmax: %f, dr: %f", rmax, dr)
    rdf = freud.density.RDF(rmax=rmax, dr=dr)
    for index, snapshot in enumerate(sdanalysis.open_trajectory(infile)):
        if index > num_frames:
            break
        rdf.accumulate(snapshot.freud_box(), snapshot.position)

    return rdf


def _structure_factor_wave_number(
    rdf: freud.density.RDF, wave_number: float, num_particles: int
):
    """Calculate the static structure factor for a single wave-number."""
    dr = rdf.R[1] - rdf.R[0]
    integral = dr * np.sum((rdf.RDF - 1) * rdf.R * np.sin(wave_number * rdf.R))
    density = num_particles / rdf.box.volume
    return 1 + 4 * np.pi * density / wave_number * integral


def static_structure_factor(infile, num_frames):
    num_particles = next(sdanalysis.open_trajectory(infile)).num_mols

    rdf = radial_distribution(infile, num_frames=num_frames)

    ssf = []
    xvalues = np.zeros_like(rdf.R)
    xvalues[:] = rdf.R
    xvalues = xvalues[2:]
    for value in xvalues:
        ssf.append(_structure_factor_wave_number(rdf, value, num_particles))

    logger.debug("x values: %s", rdf.R)
    return xvalues, ssf


@click.group()
def main():
    pass


@main.command()
@click.argument("infile", type=click.Path())
@click.argument("outfile", type=click.Path())
@click.option("--num-frames", default=10000, type=int)
def plot_rdf(infile, outfile, num_frames):
    """Plot the radial distribution function of an input trajectory."""
    rdf = radial_distribution(infile, num_frames=num_frames)

    fig, ax = plt.subplots()
    ax.plot(rdf.R, rdf.RDF)
    ax.set_xlabel(r"$r$")
    ax.set_ylabel("RDF")
    fig.savefig(outfile)


@main.command()
@click.argument("infile", type=click.Path())
@click.argument("outfile", type=click.Path())
@click.option("--num-frames", default=1000, type=int)
def plot_ssf(infile, outfile, num_frames):
    """Plot the static structure factor of an input trajectory."""
    x, ssf = static_structure_factor(infile, num_frames)
    logger.debug("x values: %s, y values %s", x, ssf)

    fig, ax = plt.subplots()
    ax.plot(x, ssf)
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$S(k)$")
    fig.savefig(outfile)


if __name__ == "__main__":
    main()
