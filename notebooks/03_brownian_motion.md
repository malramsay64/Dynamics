---
jupyter:
  jupytext:
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
import numba

from sdanalysis import dynamics

import sys

sys.path.append("../src")

from brownian import brownian
import figures
figures.use_my_theme()

alt.data_transformers.enable("csv")
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
def get_relax(steps: int = 1000, time: int = 10, step_size: float = 0.25, molecules: int = 2000):
    all_tau_L = np.zeros(molecules)
    all_tau_F = np.zeros(molecules)
    delta = np.linalg.norm(
        brownian(np.zeros((2, molecules)), steps, time / steps, step_size), axis=0
    )
    tau_F = dynamics.MolecularRelaxation(molecules, 0.4)
    tau_L = dynamics.LastMolecularRelaxation(molecules, 0.4, 1.0)
    for time in range(delta.shape[-1]):
        tau_F.add(time, delta[:, time])
        tau_L.add(time, delta[:, time])
    return tau_F, tau_L
```

```python
num_samples = 10_000
```

```python
relax = get_relax(steps=10000, time=10, step_size=0.75, molecules=num_samples)
```

```python
brownian_relation = pandas.DataFrame({"tau_F": relax[0].get_status(), "tau_L": relax[1].get_status()})
brownian_relation = brownian_relation[brownian_relation.tau_L != 2 ** 32 - 1]
```

```python
c = (
    alt.Chart(brownian_relation.sample(num_samples))
    .mark_point()
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_L", scale=alt.Scale(type="log")),
    )
)
c
```

```python
df = brownian_relation[["tau_L", "tau_F"]]
for x1, x2 in itertools.combinations(df.columns, 2):
    correlation, pValue = scipy.stats.pearsonr(getattr(df, x1), getattr(df, x2))
    print(f"{x1: <8} {x2: <8} {correlation:.2f}")
```

```python
# brownian_relation = brownian_relation.loc[:, ["tau_L", "tau_F"]]
brownian_relation["dataset"] = "Brownian"

df = pandas.read_hdf("../data/analysis/dynamics.h5", "molecular_relaxations")
df = df[df.tau_L != 2 ** 32 - 1]
df = df[df.tau_F != 2 ** 32 - 1]

lowT_mask = (df["pressure"] == 13.50) & (df["temperature"] == 1.30)
lowT_df = df[lowT_mask].copy()
lowT_df["dataset"] = "Low T"
lowT_df.reset_index(inplace=True)

highT_mask = (df["pressure"] == 13.50) & (df["temperature"] == 2.50)
highT_df = df[highT_mask].copy()
highT_df["dataset"] = "High T"
highT_df.reset_index(inplace=True)
```

```python
c = (
    alt.Chart(lowT_df.sample(num_samples))
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_L", scale=alt.Scale(type="log")),
    )
)
c
```

```python
c = (
    alt.Chart(highT_df.sample(num_samples))
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_L", scale=alt.Scale(type="log")),
    )
)
c
```

```python
df = pandas.concat(
    [
        brownian_relation.sample(num_samples),
        lowT_df.sample(num_samples),
        highT_df.sample(num_samples),
    ],
    sort=True,
)
df.reset_index(drop=True, inplace=True)
```

```python
df["diffs"] = df.tau_L - df.tau_F
df["log_diffs"] = np.log10(df.diffs[df.diffs > 0])
df["log_tau_F"] = np.log10(df.tau_F)
df=df.query("log_diffs > 0")
```

```python
alt.Chart(df).mark_bar(opacity=0.8).encode(
    x=alt.X("log_diffs:Q", bin=alt.Bin(maxbins=100)),
    y=alt.Y("count()", stack=None),
    color="dataset",
)
```

```python
sim_df = pandas.read_hdf("../data/analysis/dynamics.h5", "molecular_relaxations")
sim_df = sim_df[sim_df["pressure"] == 13.50]
sim_df = sim_df[sim_df.tau_L != 2 ** 32 - 1]
sim_df = sim_df[sim_df.tau_F != 2 ** 32 - 1]
sim_df = sim_df.reset_index()
sim_df = sim_df[["tau_F", "tau_L", "temperature"]]
```

```python
br_df = brownian_relation[["tau_F", "tau_L"]]
br_df["temperature"] = "Brownian"
```

```python
sim_df["temperature"] = sim_df["temperature"].astype(str)
```

```python
np
```

```python
all_df = pandas.concat([br_df, sim_df])
```

```python
all_df["diffs"] = all_df.tau_L - all_df.tau_F
all_df["log_diffs"] = np.nan
pos_diffs = all_df["diffs"] > 0
all_df.loc[pos_diffs, "log_diffs"] = np.log10(all_df.loc[pos_diffs, "diffs"])
all_df["log_tau_F"] = np.log10(all_df.tau_F)
all_df=all_df.query("log_diffs > 0")
```

```python
all_groups = []
for index, group in all_df.groupby("temperature"):
    all_groups.append(group.sample(num_samples, replace=True))
all_plot_df = pandas.concat(all_groups)
```

```python
alt.Chart(all_plot_df).mark_bar(opacity=0.8).encode(
    x=alt.X("log_diffs:Q", bin=alt.Bin(maxbins=100)),
    y=alt.Y("count()", stack=None),
    color="temperature",
)
```

```python

```
