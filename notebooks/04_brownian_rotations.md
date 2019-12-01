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

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas
import altair as alt
import itertools
import scipy.stats

from sdanalysis import dynamics
from dynamics_analysis import brownian, figures

```


```python
def get_rotational_relax(
    steps: int, time: int = 10, step_size: float = 0.005, molecules: int = 2000
):
    delta_theta = np.angle(
        brownian.brownian_rotation(
            np.ones(molecules, np.complex128), steps, steps / time, step_size
        )
    )
    tau_T2 = dynamics.MolecularRelaxation(molecules, np.pi / 2)
    tau_T3 = dynamics.MolecularRelaxation(molecules, np.pi / 3)
    tau_T4 = dynamics.MolecularRelaxation(molecules, np.pi / 4)
    for time in range(delta_theta.shape[-1]):
        tau_T2.add(time, delta_theta[:, time])
        tau_T3.add(time, delta_theta[:, time])
        tau_T4.add(time, delta_theta[:, time])
    return tau_T2, tau_T2, tau_T4
```


```python
num_samples = 10_000
relax = get_rotational_relax(10_000, time=10, molecules=num_samples)
df_brownian = pandas.DataFrame(
    {
        "tau_T2": relax[0].get_status(),
        "tau_T3": relax[1].get_status(),
        "tau_T4": relax[2].get_status(),
        "temperature": "Brownian",
    }
)
```

```python
alt.Chart(df_brownian).mark_bar().encode(
    x=alt.X("tau_T2", bin=alt.Bin(maxbins=100)), y="count()"
)
```

```python
df_brownian["tau_T2"].mean() / df_brownian["tau_T4"].mean()
```

```python
df = pandas.read_hdf("../data/analysis/dynamics_clean.h5", "molecular_relaxations")
df = df.query("pressure==13.50")
df = df.mask(df > 2 ** 20).dropna()
df["temperature"] = df["temperature"].astype(str)
```

```python
num_samples = 1000
df_all = pandas.concat([df, df_brownian], sort=False)
all_groups = []
for index, group in df_all.groupby("temperature"):
    group = group.sample(num_samples, replace=True).set_index("temperature")
    group = group / group["tau_T2"].mean()
    all_groups.append(group)
df_plot = pandas.concat(all_groups).reset_index()
```

```python
alt.Chart(df_plot).mark_line(opacity=0.8).encode(
    x=alt.X("tau_T4", title="tau_T2 / <tau_T2>", bin=alt.Bin(maxbins=50)),
    y=alt.Y("count()", stack=None),
    color="temperature",
)
```

```python
np.mean(df["tau_T2"]) / scipy.stats.hmean(df["tau_T2"])
```

```python
df_brownian["tau_T2"].agg(["mean", lambda x: scipy.stats.hmean(x.values)])
```

```python
df_plot.groupby("temperature")["tau_T2"].agg(
    [lambda x: np.mean(x) / scipy.stats.hmean(x.values)]
)
```

```python
alt.Chart(df_plot).mark_line(opacity=0.8).transform_calculate(
    "ratio", alt.datum.tau_T2 / alt.datum.tau_T4
).encode(
    x=alt.X("ratio:Q", title="tau_T2 / <tau_T2>", bin=alt.Bin(maxbins=20)),
    y=alt.Y("count()", stack=None),
    color="temperature",
)
```

```python
df_plot
```
