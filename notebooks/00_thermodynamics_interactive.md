---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.2.4
  kernelspec:
    display_name: dynamics
    language: python
    name: dynamics
---

# Visualisation of Thermodynamic Quantities

One of the first steps for determining a simulation has proceeded as expected,
particularly those that are supposed to be at equilibrium,
is to check the thermodynamic quantities.
For simulations at equilibrium,
the thermodynamics should remain relatively constant
over the course of the simulation.

This provides an interactive method of investigating
all the values to check that they are as expected.

```python
from pathlib import Path

import pandas as pd
import ipywidgets as widgets
import matplotlib.pyplot as plt
import altair as alt

from sdanalysis.util import get_filename_vars
from sdanalysis.figures.interactive_config import parse_directory

import sys

sys.path.append("../src")
import figures
```

```python
dataset = parse_directory(
    directory=Path("../data/simulations/trimer/output/"), glob="thermo-*log"
)
```

```python
discard_quantities = ["timestep", "lz", "xz", "yz", "N"]
file = "../data/simulations/trimer/output/thermo-Trimer-P13.50-T1.25.log"
with open(file) as src:
    thermo_quantities = [c.strip() for c in src.readline().split("\t")]
    thermo_quantities = [c for c in thermo_quantities if c not in discard_quantities]
    thermo_quantities = [c for c in thermo_quantities if "rigid_center" not in c]
```

```python
pressure = widgets.ToggleButtons(description="Pressure", options=list(dataset.keys()))
temperature = widgets.ToggleButtons(
    description="Temperature", options=list(dataset.get(pressure.value).keys())
)

quantities = widgets.ToggleButtons(description="Quantity", options=thermo_quantities)


def update_temperatures(change):
    temperature.options = list(dataset.get(change.new).keys())


pressure.observe(update_temperatures, names="value")


@widgets.interact(pressure=pressure, temperature=temperature, quantity=quantities)
def plot_figure(pressure, temperature, quantity):
    filename = dataset.get(pressure).get(temperature).get("None").get("None")
    df = pd.read_csv(
        filename,
        sep="\t",
        index_col="timestep",
        usecols=["timestep", quantity],
        converters={"timestep": pd.to_timedelta},
    )
    df = df[~df.index.duplicated(keep="first")]
    df = df.resample("10ms").agg(["mean", "std"])
    df.columns = [col[-1] for col in df.columns.values]
    df = df.reset_index()
    df["time"] = df["timestep"].astype(int) * 0.005
    df = df.drop(columns="timestep")
    c = alt.Chart(df).encode(
        x=alt.X("time", title="Time", axis=alt.Axis(format="e")),
        y=alt.Y("mean", title=quantity, scale=alt.Scale(zero=False)),
    )
    c = c.mark_line() + c.mark_errorband().encode(yError=alt.YError("std"))
    return c
```

```python

```
