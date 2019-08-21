#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

"""Code to assist in the creation of figures."""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import altair as alt
import click
import freud
import matplotlib.pyplot as plt
import numpy as np
import pandas
import sdanalysis

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def my_theme() -> Dict[str, Any]:
    """Define an Altair theme to use in all visualisations.

    This defines a simple theme, specifying the aspect ratio of 4:6 
    and removing the grid from the figure which is distracting.

    """
    return {
        "config": {
            "view": {"height": 400, "width": 600},
            "legend": {"titleFontSize": 20, "labelFontSize": 16},
            "axis": {"grid": False, "labelFontSize": 16, "titleFontSize": 20},
            "header": {"titleFontSize": 22, "labelFontSize": 18},
            "background": "white",
        }
    }


def json_dir(data, data_dir="altairdata"):
    data_dir = Path(data_dir)
    data_dir.mkdir(exist_ok=True)
    return alt.pipe(
        data, alt.to_json(filename=str(data_dir / "{prefix}-{hash}.{extension}"))
    )


def use_my_theme():
    """Register and my custom Altair theme."""
    # register and enable the theme
    alt.themes.register("my_theme", my_theme)
    alt.themes.enable("my_theme")


def use_data_transformer():
    """Register and use an altair data transformer"""
    alt.data_transformers.register("json_dir", json_dir)
    alt.data_transformers.enable("json_dir")


# This configures Altair to behave as expected when I load this module
use_my_theme()
use_data_transformer()


def plot_dynamics(
    df: pandas.DataFrame, prop: str, title: Optional[str] = None, scale: str = "linear"
) -> alt.Chart:
    """Helper to plot dynamics quantities using Altair.

    Args:
        df: DataFrame containing quantities to plot. There should be columns 
            with suffixes '_value', '_lower', and '_upper'.
        prop: The property to plot, with the respective values for that property 
            having suffixes.
        title: Custom title used to label the property.
        scale: Type of scale to use for the property, this should be "linear" or "log".

    """
    if title is None:
        title = prop
    axis_format = "g"
    if scale == "log":
        axis_format = "e"

    dyn_chart_base = (
        alt.Chart(df)
        .encode(
            x=alt.X(
                "time",
                title="Timesteps",
                scale=alt.Scale(type="log"),
                axis=alt.Axis(format="e"),
            ),
            color=alt.Color("temperature:N", title="Temperature"),
        )
        .transform_filter(alt.datum.msd_value < 50)
    )

    confidence_interval = dyn_chart_base.mark_errorband().encode(
        y=alt.Y(
            prop + "_lower:Q",
            title=title,
            scale=alt.Scale(type=scale),
            axis=alt.Axis(format=axis_format),
        ),
        y2=alt.Y2(prop + "_upper:Q", title=title),
    )

    values = dyn_chart_base.mark_line().encode(
        y=alt.Y(
            prop + "_value:Q",
            title=title,
            scale=alt.Scale(type=scale),
            axis=alt.Axis(format=axis_format),
        )
    )

    return values + confidence_interval


def plot_relaxations(
    df: pandas.DataFrame, prop: str, title: Optional[str] = None
) -> alt.Chart:
    """Helper to plot relaxation quantities using Altair.

    Args:
        df: DataFrame containing quantities to plot. There should be columns 
            with suffixes '_value', '_lower', and '_upper'.
        prop: The property to plot, with the respective values for that property 
            having suffixes.
        title: Custom title used to label the property.

    """
    if title is None:
        title = prop
    axis_format = "e"

    relax_chart_base = alt.Chart(df).encode(
        x=alt.X("temp_norm:Q", title="Tₘ/T", axis=alt.Axis(format="g")),
        color=alt.Color("pressure:N", title="Pressure"),
    )

    confidence_interval = relax_chart_base.mark_rule(opacity=1.0).encode(
        y=alt.Y(
            prop + "_lower:Q",
            title=title,
            scale=alt.Scale(type="log"),
            axis=alt.Axis(format=axis_format),
        ),
        y2=alt.Y2(prop + "_upper:Q", title=title),
    )

    values = relax_chart_base.mark_point().encode(
        y=alt.Y(
            prop + "_value:Q",
            title=title,
            scale=alt.Scale(type="log"),
            axis=alt.Axis(format=axis_format),
        )
    )

    return values + confidence_interval


def reshape_dataframe(df: pandas.DataFrame) -> pandas.DataFrame:
    values = []
    columns = []
    r_df = df.set_index(["temperature", "pressure", "temp_norm"])
    for col_name in r_df.columns:
        col_split = col_name.split("_")
        columns.append("_".join(col_split[:-1]))
        values.append(col_split[-1])
    r_df.columns = pandas.MultiIndex.from_arrays([columns, values])
    return (
        r_df.stack(level=0)
        .reset_index()
        .rename(index=str, columns={"level_3": "variable"})
    )


def plot_multi_relaxations(
    df: pandas.DataFrame, prop: List[str], title: str = "Relaxation Values"
) -> alt.Chart:
    """Helper to plot relaxation quantities using Altair.

    Args:
        df: DataFrame containing quantities to plot. There should be columns 
            with suffixes '_value', '_lower', and '_upper'.
        prop: The property to plot, with the respective values for that property 
            having suffixes.
        title: Custom title used to label the property.

    """
    if isinstance(prop, str):
        prop = [prop]

    axis_format = "e"

    #  plot_df = df.query("variable in @prop")
    #  plot_df.loc[:, "variable"] = title

    relax_chart_base = (
        alt.Chart(df)
        .encode(
            x=alt.X("temp_norm:Q", title="Tₘ/T", axis=alt.Axis(format="g")),
            color=alt.Color("pressure:N", title="Pressure"),
        )
        .transform_filter(alt.FieldOneOfPredicate(field="variable", oneOf=prop))
    )

    confidence_interval = relax_chart_base.mark_rule(opacity=1.0).encode(
        y=alt.Y(
            "lower:Q",
            title=title,
            scale=alt.Scale(type="log"),
            axis=alt.Axis(format=axis_format),
        ),
        y2=alt.Y2("upper:Q", title=title),
    )

    values = relax_chart_base.mark_point().encode(
        y=alt.Y(
            "value:Q",
            title=title,
            scale=alt.Scale(type="log"),
            axis=alt.Axis(format=axis_format),
        ),
        shape=alt.Shape("variable", title="Relaxation"),
    )

    return values + confidence_interval


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
