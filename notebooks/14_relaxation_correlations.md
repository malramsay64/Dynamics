---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.2.1
  kernelspec:
    display_name: dynamics
    language: python
    name: dynamics
---

# Correlation of Relaxation Times

This is a notebook to understand how relaxation times correlate with each other,
in particular investigating the relaxation of individual molecules.

```python
from pathlib import Path
import numpy as np
import scipy.stats
import pandas
from bokeh.plotting import figure, show, output_notebook
from bokeh.palettes import Dark2
from bokeh.transform import factor_cmap
from bokeh.palettes import Spectral6
import ipywidgets
import itertools
from IPython.display import display
import altair as alt

import sys
sys.path.append("../src")
import figures

figures.use_my_theme()

output_notebook()

alt.data_transformers.enable("json")
```

## Data Import

The data used for this analysis is from a series of simulations studying the dynamics of a Trimer molecule at a pressure of 13.50.
These simulations are all of equilibrated systems at a range of temperatures.
For the purposes of this investigation only a single temperature has been studied,
although it is feasible to extend this to multiple temperatures.

```python
# Modify these variables to change the simulation being investigated
temperature = "1.35"
pressure = "13.50"

# Select directory data stored and the file within that directory
infile = "../data/analysis/dynamics.h5"

# Load input file into a pandas dataframe
df = pandas.read_hdf(infile, "molecular_relaxations")

# There are a number of initial configurations in the data, only select the first
# df = pandas.concat(chunks)
```

## Visualisation

A simple method of understanding how these timescales relate to one another
is to plot them against each other in a scatter plot.
Where there is some correlation between them this will exist on a straight line throughout the figure.
Since the values of the timescales I am investigating are over such a large range,
both the x and y axis have a log scale.

The timescales I am plotting in this figure are the time for a molecule to undergo a certain motion,
with all motions relative to the initial position.
The motions are tabulated below;

- `tau_F` -> travel a distance of 0.4
- `tau_D`  -> travel a distance of 1
- `tau_L` -> last passage of a distance of 0.4 before reaching 1.0
- `tau_T2`  -> rotate a distance of $\pi/2$
- `tau_T3`  -> rotate a distance of $\pi/3$
- `tau_T4`  -> rotate a distance of $\pi/4$

```python
# Create buttons to select the quantity plotted on each axis
x_opts = ipywidgets.ToggleButtons(options=df.columns, description="X axis:")
y_opts = ipywidgets.ToggleButtons(
    value=df.columns[1], options=df.columns, description="Y axis:"
)


@ipywidgets.interact(x=x_opts, y=y_opts)
def create_plot(x, y):
    f = figure(x_axis_type="log", y_axis_type="log", x_axis_label=x, y_axis_label=y)
    f.scatter(x, y, source=df)
    show(f)
```

## Quantitation of correlations

While it is nice to look at a picture,
it is difficult to make a fair comparison between two images.
Here I aim to quantify the level of correlation between the timescales.
The metric I have chosen is the [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
which is implemented in the `scipy.stats` module as [pearsonr](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.pearsonr.html).

```python
for x1, x2 in itertools.combinations(df.columns, 2):
    correlation, pValue = scipy.stats.pearsonr(getattr(df, x1), getattr(df, x2))
    print(f"{x1: <8} {x2: <8} {correlation:.2f}")
```

It is notable that the correlation between the last passage time `tau_L` and the diffusion time `tau_D` are so similar.
This really is indicating that the last passage time is a good proxy for diffusive behaviour.
There is still a strong correlation between the last passage time and the relaxation times,
so there is a reasonable argument to be made that these values are all connected in some way.

```python
selected = ["tau_L", "tau_T2", "tau_T3", "tau_T4"]
data = df[selected].melt("tau_L")
f = figure(
    x_axis_type="log", y_axis_type="log", x_axis_label="tau_L", y_axis_label="tau_Tn"
)
f.scatter(
    "tau_L",
    "value",
    source=data,
    color=factor_cmap(
        "variable", palette=Spectral6, factors=["tau_T2", "tau_T3", "tau_T4"]
    ),
)
show(f)
```

```python
f = figure(
    x_axis_type="linear",
    y_axis_type="linear",
    x_axis_label="tau_L",
    y_axis_label="tau_Tn",
)
for col, colour in zip(["tau_T2", "tau_T3", "tau_To"], Dark2[3]):
    f.scatter(
        "tau_Data", col, color=colour, legend=col[4:], source=df, alpha=0.6, size=3
    )
f.legend.location = "top_left"
f.legend.click_policy = "hide"
show(f)
```

```python
data_log = data.copy()
data_log.value = np.log10(data_log.value / data_log["tau_L"])
alt.Chart(data_log).mark_area(opacity=0.7, interpolate="step").encode(
    alt.X("value", bin=alt.Bin(maxbins=100)),
    alt.Y("count()", stack=None),
    alt.Color("variable"),
)
```

```python

```
