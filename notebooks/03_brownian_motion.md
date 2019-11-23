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

# Brownian Motion

This is an attempt at a better understanding of our new metrics of molecular motion,
namely how they relate to each other when Brownian Dynamics are observed.
This notebook implements Brownian dynamics using the recipe from the [scipy cookbook]( http://scipy-cookbook.readthedocs.io/items/BrownianMotion.html),
then uses the simulation of Brownian motion to investigate
how the molecular relaxation times respond.


## Implementation

The code in the cell below implements the Brownian dynamics.
For 2D Brownian dynamics, x0 with 2 elements can be used as the input.

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas
import altair as alt
import itertools
import scipy.stats

from sdanalysis import dynamics
from dynamics_analysis import figures, brownian

```

### Testing Brownian Dynamics

This is a simple test case of the function to Generate the Brownian dynamics
ensuring that the resulting trajectory is sensible,
and to give some idea of the spread of values.

```python
x = brownian.brownian(np.zeros((2,)), 100, 0.25, 0.25)
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

## Dynamics Analysis

Now we have a function to generate Brownian dynamics,
I want to use it to understand the behaviour
of the new molecular quantities I have defined.
The calculation of these relaxation quantities
is done using the `get_relax` function
from the `brownian` module.

```python
num_samples = 20_000
relax = brownian.get_brownian_relax(
    steps=10_000, time=10, step_size=0.75, molecules=num_samples
)
```

Once these relaxation values are calculated,
there is some clean-up that needs to occur.
Firstly the conversion to a `pandas.DataFrame`
which allows for the easy manipulation of the data.
In this form,
the missing values need to be removed.
Any values which did not complete relaxation
have a value of $2^{32} - 1$,
the largest value which can be stored as a `uint32`
so these values are removed.
The final pre-processing step is a timescale normalisation.
To be able to compare the Brownian motion
with those from simulations
I am normalising by the mean value of
the structural relaxation time $\tau_F$.

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

One of the results of Brownian motion
is that there should be no correlation
between particles on short and long timescales,
that it particles with a fast structural relaxation (short time)
shouldn't correlate with particles that have a fast diffusion (long time).

```python
c_brownian = (
    alt.Chart(df_brownian.sample(num_samples, replace=True))
    .mark_point(opacity=0.3)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_D", scale=alt.Scale(type="log")),
    )
)
with alt.data_transformers.enable("default"):
    c_brownian.save("figures/brownian_tau_F_tau_D.svg", webdriver="firefox")
```

This relation does indeed hold true here as shown in the below figure.

![tau_F vs tau_D](figures/brownian_tau_F_tau_D.svg)

One method we can use to measure the correlation
of these two values is through the Pearson Correlation Coefficient,
which is a measure of how well the ordering of one collection of values
corresponds to another.
The Pearson correlation has values between -1 and 1,
with a value of -1 indicating anti-correlation,
that is the fastest structural relaxation is the slowest diffuser,
a value of 0 indicating no correlation,
and a value of 1 indicating perfect correlation.

```python
correlation, pValue = scipy.stats.pearsonr(df_brownian["tau_F"], df_brownian["tau_D"])
print(f"Correlation of tau_F and tau_D {correlation:.2f}")
```

Using the scipy implementation of the Pearson correlation,
we get a value very close to 0,
meaning there is no correlation between these two quantities
as expected.

## Comparison to Molecular System

Using these tools to analyse brownian motion
is a great test of their applicability,
however it doesn't provide us with any additional information.
For that we need to study a molecular system.

The data for the molecular systems is within the file
`../data/analysis/dynamics_clean.h5` under the `molecular_relaxations` table.
This data is calculated when running the command `make dynamics`.

There are two separate datasets which I am investigating
for the analysis of the molecular relaxations.
A high temperature simulation,
from simulations at a temperature of 2.50,
which is expected to be Brownian in nature.
Additionally there is a low temperature dataset
from a temperature of 1.30,
which is below the melting point of 1.35
and so is expected to show behaviour of a supercooled liquid.

To ensure a better comparison between the datasets
I am ensuring they all have the same number of data points
by using a method of sampling with replacement.

```python
df_brownian["dataset"] = "Brownian"

df = pandas.read_hdf("../data/analysis/dynamics_clean.h5", "molecular_relaxations")
df = df.set_index(["pressure", "temperature"]).sort_index()

df_low = df.loc[(13.50, 1.30), :].copy()
df_low = df_low / df_low["tau_F"].mean()
df_low["dataset"] = "Low T"
df_low.reset_index(inplace=True)

df_high = df.loc[(13.50, 2.50), :].copy()
df_high = df_high / df_high["tau_F"].mean()
df_high["dataset"] = "High T"
df_high.reset_index(inplace=True)

df_sim = pandas.concat(
    [
        df_low.sample(num_samples, replace=True),
        df_high.sample(num_samples, replace=True),
    ]
)
```

Plotting the data in the same way as the Brownian motion
we get very different pictures for
the low and the high temperature simulations.

```python
c = (
    alt.Chart(df_sim)
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_F", scale=alt.Scale(type="log")),
        y=alt.Y("tau_D", scale=alt.Scale(type="log")),
        column="dataset",
    )
)
with alt.data_transformers.enable("default"):
    c.save("../figures/trimer_tau_F_tau_D.svg", webdriver="firefox")
```

![tau_F vs tau_D for the trimer molecule](../figures/trimer_tau_F_tau_D.svg)

Also useful in this discussion are the correlations coefficients.

```python
correlation, pValue = scipy.stats.pearsonr(df_high["tau_F"], df_high["tau_D"])
print(f"Correlation of tau_F and tau_D {correlation:.2f}")
```

```python
correlation, pValue = scipy.stats.pearsonr(df_low["tau_F"], df_low["tau_D"])
print(f"Correlation of tau_F and tau_D {correlation:.2f}")
```

The value of tau_F for the low temperature dataset
is distributed over a much wider range of values
than of the high temperature dataset.
While there is roughly the spread of diffusive motion.
The distribution of tau_F is representative of
the dynamic heterogeneities present within
the low temperature simulation.
It is important to note
that the dynamic heterogeneities
are only present on the short
time and length scales.

## Reversals

One of the results of comparing dynamics with the standard quantities
is that the first relaxation is not indicative of relaxation,
instead the last passage time is more important.
This can be explained by the reversal of structural relaxations
which are captured by the traditional quantities.

The distances of reversals which we are discussing here
being 0.43 are beyond what would typically be considered reversible.
That is, the molecule has escaped the confinement of it's local area.
It is important to note that there is also
a rate at which reversals occur as modelled by Brownian motion.
In fact when considering a random walk,
in 2D it is guaranteed you will pass back through the origin.

When modelling the last passage time using Brownian motion.

```python
c = (
    alt.Chart(df_brownian)
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_L", scale=alt.Scale(type="log")),
        y=alt.Y("tau_D", scale=alt.Scale(type="log")),
    )
)
with alt.data_transformers.enable("default"):
    c.save("../figures/brownian_tau_L_tau_D.svg", webdriver="firefox")
```

![tau_L vs tau_D for brownian relaxation](../figures/brownian_tau_L_tau_D.svg)

Which translates to the following in simulations

```python
c = (
    alt.Chart(df_sim)
    .mark_point(opacity=0.2)
    .encode(
        x=alt.X("tau_L", scale=alt.Scale(type="log")),
        y=alt.Y("tau_D", scale=alt.Scale(type="log")),
        column="dataset",
    )
)
with alt.data_transformers.enable("default"):
    c.save("../figures/trimer_tau_L_tau_D.svg", webdriver="firefox")
```

![tau_L vs tau_D for the trimer molecule](../figures/trimer_tau_L_tau_D.svg)

Also useful in this discussion are the correlations coefficients.

```python
correlation, pValue = scipy.stats.pearsonr(df_high["tau_L"], df_high["tau_D"])
print(f"Correlation of tau_L and tau_D {correlation:.2f}")
```

```python
correlation, pValue = scipy.stats.pearsonr(df_low["tau_L"], df_low["tau_D"])
print(f"Correlation of tau_L and tau_D {correlation:.2f}")
```

There is a strong correlation between these two quantities,
an element of which is to be expected,
since the diffusion length is included as part of the last passage time,
defining the last part.
Despite this, it is interesting to compare with the first passage time.
This is directly a measure of the idea of cage escape,
for the molecule to undergo diffusive motion it has to have escaped the cage,
1.2 units is a large motion.

Also particularly interesting here
is that the simulations have a more pronounced
relationship than the Brownian motion,
likely the effect of caging is influencing this.

### Understanding reversals

Taking the time between the first and last relaxation time
provides a nice picture.

```python
df = pandas.concat([df_brownian, df_sim], sort=True)
df["diffs"] = df.tau_D - df.tau_L
mask = df["diffs"] != 0  # Problematic with log
df.loc[mask, "log_diffs"] = np.log10(df.loc[mask, "diffs"])
df["log_tau_F"] = np.log10(df.tau_F)
# df = df.query("log_diffs > 0")
```

```python
alt.Chart(df).mark_bar(opacity=0.8).encode(
    x=alt.X("diffs:Q", bin=alt.Bin(maxbins=50)),
    y=alt.Y("count()", stack=None),
    color="dataset",
)
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
