#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Code to assist in the creation of figures."""

import altair as alt
import matplotlib.pyplot as plt
import freud
import sdanalysis
import click


def my_theme():
    return {"config": {"view": {"height": 400, "width": 600}, "axis": {"grid": False}}}


def use_my_theme():
    # register and enable the theme
    alt.themes.register("my_theme", my_theme)
    alt.themes.enable("my_theme")

@click.command()
@click.argument("infile", type=click.Path())
@click.argument("outfile", type=click.Path())
@click.option("--num-frames", default=10_000, type=int)
def plot_radial_distribution(infile, outfile, num_frames) -> plt.Figure:
    dr = 0.1
    rmax = 20
    rdf = freud.density.RDF(rmax=rmax, dr=dr)

    for index, snapshot in enumerate(sdanalysis.open_trajectory(infile)):
        if index > num_frames:
            break
        rdf.accumulate(snapshot.freud_box(), snapshot.position)

    fig, ax = plt.subplots()
    ax.plot(rdf.R, rdf.RDF)
    ax.set_xlabel("r")
    ax.set_ylabel("RDF")

    fig.savefig(outfile)

    return fig

if __name__ == "__main__":
    plot_radial_distribution()
