---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.3.0
  kernelspec:
    display_name: dynamics
    language: python
    name: dynamics
---

# Computing Confidence Intervals

The computing of confidence intervals is an important part of a simulation,
it provides information on the level of accuracy of a particular analysis.
Since most of the distributions of values present in the dynamics
are not normally distributed,
that is they are far from that of a standard normal distribution,
instead of using standard statistical tools
for which the theory deriving them assumes a normal distribution,
I will use the bootstrapping method.
Bootstrapping is a Monte Carlo method of calculating statistics,
randomly sampling form the measured distribution
to generate a mean value and a confidence interval.

```python
# Make dealing with filesystem paths much, much simpler
from pathlib import Path

# Read/write data files and data analysis
import pandas
import numpy as np
import altair as alt
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats

from sdanalysis.relaxation import series_relaxation_value

from dynamics_analysis import figures

```

For this analysis we are going to be using
the data from the cleaned dynamics file.
This data has all the values collated from
the individual temperature and pressure analyses
and has then the values for which there are
only a small number of observations have been removed.

```python
# Where the data files with the results are located
data_dir = Path("../data/analysis")

# Load data for most of the figures
dynamics_df = pandas.read_hdf(data_dir / "dynamics_clean.h5", "dynamics")
```

To understand the bootstrapping procedure
we want just a single observation,
so for this we group the dataset into
the component temperatures, pressures, and time intervals.
This breakup leaves each group with many
independent observations from a range of stating configurations.
It is these independent observations
we can use for the bootstrapping.

```python
index, group = next(iter(dynamics_df.groupby(["temperature", "pressure", "time"])))
group.head()
```

The reason we are using the bootstrap method
is that the measured values
are not normally distributed.
We can see this in the mean squared displacement (MSD)
plotted below.
Instead of a normal distribution,
where values are evenly distributed on either side
instead we have a highly skewed distribution of values.

```python
alt.Chart(group).mark_bar().encode(alt.X("rot2", bin=alt.Bin(maxbins=50)), y="count()")
```

The bootstrapping procedure generates
a distribution of results
in a method demonstrated below.
The measured distribution of values `group.rot2`
is sampled the with replacement
for the number of observations made.
This is repeated many times (in this case 1000)
taking the mean of each of the 1000 samples.

```python
results = np.random.choice(
    # The distribution of values to sample
    group.rot2,
    # Shape of the sampled values,
    # being the length of the distribution x 1000 replications
    (len(group), 1000),
    # The sampling is performed with replacement
    replace=True,
).mean(
    axis=1
)  # Find the mean value for each replication
```

This creates a distribution of values
from which we can extract a confidence interval
from the range of possible values
and a mean from the median of the results.

```python
alt.Chart(pandas.DataFrame({"results": results})).mark_bar().encode(
    alt.X("results", bin=alt.Bin(maxbins=20)), y="count()"
)
```

This procedure is packaged up in the bootstrapped package,
and the mean with confidence interval can be calculated
as shown below.

```python
val = bs.bootstrap(group.rot2.values, bs_stats.mean)
val
```

## Confidence Interval of Curve Fits

Finding the parameters for a curve fit
through bootstrapping is a similar procedure.
However in this case the independent observations
are the parameters of the curve fit
for each of the starting configurations.

```python
val_index, df = next(iter(dynamics_df.groupby(["temperature", "pressure"])))
df = df.set_index("time")
index, group = next(iter(df.groupby("keyframe")))
group.head()
```

```python
series_relaxation_value(group.rot2)
```

```python
from sdanalysis.relaxation import series_relaxation_value
```

```python
results = df.groupby("keyframe").struct.agg(series_relaxation_value)
```

```python
alt.Chart(pandas.DataFrame({"results": results})).mark_bar().encode(
    alt.X("results", bin=alt.Bin(maxbins=20)), y="count()"
)
```

```python
df_relax = df.groupby("keyframe").agg(series_relaxation_value)
val = bs.bootstrap(df_relax.struct.values, bs_stats.mean)
val
```

```python
df_relax = dynamics_df.set_index("time")
df_agg = df_relax.groupby(["temperature", "pressure", "keyframe"]).agg(
    series_relaxation_value
)
```

```python
df_vals = df_agg.groupby(["temperature", "pressure"]).agg("mean")
```

```python
alt.Chart(df_vals.reset_index().query("pressure == 13.50")).mark_point().encode(
    x="temperature", y="rot2"
)
```
