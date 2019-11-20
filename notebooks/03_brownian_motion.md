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

Brownian Motion
===============

This is an attempt at a better understanding of our new metrics of molecular motion,
namely how they relate to each other when Brownian Dynamics are observed.
This notebook implements Brownian dynamics using the recipe from the [scipy cookbook]( http://scipy-cookbook.readthedocs.io/items/BrownianMotion.html),
then uses the simulation of Brownian motion to investigate
how the molecular relaxation times respond.


Implementation
--------------

The code in the cell below implements the brownian dynamics.
For 2D brownian dynamics, x0 with 2 elements can be used as the input.

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas
import altair as alt
import itertools
import scipy.stats

from sdanalysis import dynamics

import sys

sys.path.append("../src")

from brownian import brownian
import figures

```

### Testing Brownian Dynamics

This is a simple test case of the function to Generate the Brownian dynamics
ensuring that the resulting trajectory is sensible,
and to give some idea of the spread of values.

```python
x = brownian(np.zeros((2,)), 100, 0.25, 0.25)
```

```python
# Plot the 2D trajectory.
plt.plot(x[0], x[1])

# Mark the start and end points.
plt.plot(x[0, 0], x[1, 0], "go")
plt.plot(x[0, -1], x[1, -1], "ro")

# More plot decorations.
plt.title("2D Brownian Motion")
plt.xlabel("x", fontsize=16)
plt.ylabel("y", fontsize=16)
plt.axis("equal")
plt.grid(True)
plt.show()
```

Dynamics Analysis
-----------------

Now we have a function to generate brownian dynamics,
I want to use it to




```python
def get_relax(
    steps: int = 1000, time: int = 10, step_size: float = 0.25, molecules: int = 2000
):
    delta = np.linalg.norm(
        brownian(np.zeros((2, molecules)), steps, time / steps, step_size), axis=0
    )
    tau_F = dynamics.MolecularRelaxation(molecules, 0.4)
    tau_D = dynamics.MolecularRelaxation(molecules, 2.0)
    tau_L = dynamics.LastMolecularRelaxation(molecules, 0.4, 2.0)
    for time in range(delta.shape[-1]):
        tau_F.add(time, delta[:, time])
        tau_D.add(time, delta[:, time])
        tau_L.add(time, delta[:, time])
    return tau_F, tau_D, tau_L
```


```python
num_samples = 20_000
relax = get_relax(steps=10000, time=10, step_size=0.75, molecules=num_samples)
```

```python
df_brownian = pandas.DataFrame(
    {
        "tau_F": relax[0].get_status(),
        "tau_D": relax[1].get_status(),
        "tau_L": relax[2].get_status(),
    }
)
df_brownian = df_brownian.mask(df_brownian == 2 ** 32 - 1).dropna()
df_brownian = df_brownian / df_brownian["tau_F"].mean()
```

```python
c_brownian = (
    alt.Chart(df_brownian.sample(num_samples, replace=True))
    .mark_point(opacity=0.3)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_D", scale=alt.Scale(type="log")),
    )
)
c_brownian
```

```python
df = df_brownian[["tau_D", "tau_L"]]
for x1, x2 in itertools.combinations(df.columns, 2):
    correlation, pValue = scipy.stats.pearsonr(getattr(df, x1), getattr(df, x2))
    print(f"{x1: <8} {x2: <8} {correlation:.2f}")
```

```python
# brownian_relation = brownian_relation.loc[:, ["tau_L", "tau_F"]]
df_brownian["dataset"] = "Brownian"

df = pandas.read_hdf("../data/analysis/dynamics_clean.h5", "molecular_relaxations")
df = df.set_index(["pressure", "temperature"]).sort_index()

df_low = df.loc[(13.50, 1.30), :].copy()
df_low = df_low / df_low["tau_F"].mean()
df_low["dataset"] = "Low T"
df_low.reset_index(inplace=True)

df_high = df.loc[(13.50, 2.50), :].copy()
df_high = df_high.query("tau_F < 1e7")
df_high = df_high / df_high["tau_F"].mean()
df_high["dataset"] = "High T"
df_high.reset_index(inplace=True)
```

```python
c_low = (
    alt.Chart(df_low.sample(num_samples, replace=True))
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_D", scale=alt.Scale(type="log")),
    )
)
c_low
```

```python
c_high = (
    alt.Chart(df_high.sample(num_samples, replace=True))
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_L", scale=alt.Scale(type="log")),
    )
)
c_high
```

```python
df = pandas.concat(
    [
        df_brownian.sample(num_samples, replace=True),
        df_low.sample(num_samples, replace=True),
        df_high.sample(num_samples, replace=True),
    ],
    sort=True,
)
df.reset_index(drop=True, inplace=True)
```

```python
alt.Chart(df).mark_point(opacity=0.3).encode(
    x=alt.X("tau_D", scale=alt.Scale(type="log")),
    y=alt.Y("tau_L", scale=alt.Scale(type="log")),
    color="dataset",
    row="dataset",
)
```

```python
df["diffs"] = (df.tau_L - df.tau_F) / df.tau_F
df["log_diffs"] = np.log10(df.diffs[df.diffs > 0])
df["log_tau_F"] = np.log10(df.tau_F)
df = df.query("log_diffs > 0")
```

```python
alt.Chart(df).mark_bar(opacity=0.8).encode(
    x=alt.X("log_diffs:Q", bin=alt.Bin(maxbins=100)),
    y=alt.Y("count()", stack=None),
    color="dataset",
)
```

```python
df_sim = pandas.read_hdf("../data/analysis/dynamics_clean.h5", "molecular_relaxations")
df_sim = df_sim[df_sim["pressure"] == 13.50]
df_sim = df_sim.reset_index()
df_sim = df_sim[["tau_F", "tau_L", "tau_D", "temperature"]]
df_sim["temperature"] = df_sim["temperature"].astype(str)
```

```python
df_brownian["temperature"] = "Brownian"
```

```python
df = pandas.concat([df_brownian, df_sim], sort=True)
df["diffs"] = (df.tau_L - df.tau_F) / df.tau_F
df["log_diffs"] = np.nan
pos_diffs = df["diffs"] > 0
df.loc[pos_diffs, "log_diffs"] = np.log10(df.loc[pos_diffs, "diffs"])
# df_all = df_all.query("log_diffs > 0")
```

```python
all_groups = []
for index, group in df.groupby("temperature"):
    all_groups.append(group.sample(num_samples, replace=True))
df_plot = pandas.concat(all_groups)
```

```python
alt.Chart(df_plot).mark_line().encode(
    x=alt.X("log_diffs:Q", bin=alt.Bin(maxbins=50)),
    y=alt.Y("count()", stack=None),
    color="temperature",
)
```

```python
df["diffs"] = (df.tau_D - df.tau_F) / df.tau_F
df["log_diffs"] = np.nan
pos_diffs = df["diffs"] > 0
df.loc[pos_diffs, "log_diffs"] = np.log10(df.loc[pos_diffs, "diffs"])
# all_df = df.query("log_diffs > 0")

all_groups = []
for index, group in df.groupby("temperature"):
    all_groups.append(group.sample(num_samples, replace=True))
all_plot_df = pandas.concat(all_groups)

alt.Chart(all_plot_df).mark_line().encode(
    x=alt.X("log_diffs:Q", bin=alt.Bin(maxbins=50)),
    y=alt.Y("count()", stack=None),
    color="temperature",
)
```

```python
all_plot_df.groupby("temperature")["tau_D"].agg(
    [lambda x: np.mean(x) / scipy.stats.hmean(x.values)]
)
```
