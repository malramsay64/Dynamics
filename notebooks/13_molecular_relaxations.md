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

Molecular Relaxations
===========

Rather than investigating relaxations as an average over all the molecules, we want to investigate individual molecules and how their behaviour is related to that of the aggregate.

All the simulation data is from a set of simulations stored on RDS
in the folder `Dynamics/data/simulations/trimer`.
The set of simulations were all run at a pressure of 13.50
and the data analysed is collected from timesteps on an exponential scale.


## Setup

This imports the modules necessary for running the code in the rest of this notebook
while also setting up some helper functions to make the rest of the analysis simpler.

```python
import sys

import pandas as pd
import numpy as np
import altair as alt

sys.path.append("../src")
import figures

figures.use_my_theme()

alt.data_transformers.enable("json")
```

### Import Data

This reads all the data from disk,
creating a pandas DataFrame and
a Holoviews Table for analysis and plotting of the data.

```python
data_file = "../data/analysis/dynamics_clean_agg.h5"
```

```python
relaxations = pd.read_hdf(data_file, "relaxations")
relaxations = relaxations.query("pressure == 13.50")
relaxations["inv_temp"] = 1 / relaxations.temperature
```

```python
long_relax = relaxations.drop(columns=["diffusion_constant"]).melt(
    id_vars=["temperature", "pressure", "inv_temp"], var_name="Quantity"
)

# Drop invalid and abnormal values
mask = long_relax.value < 0
long_relax.value.values[mask] = np.nan
long_relax = long_relax.replace([np.inf, -np.inf], np.nan).dropna()
```

```python
relax_chart = (
    alt.Chart(long_relax)
    .mark_point(filled=True, size=100)
    .encode(
        alt.X("inv_temp", scale=alt.Scale(zero=False), axis=alt.Axis(title="1/T")),
        alt.Y(
            "value", scale=alt.Scale(type="log"), axis=alt.Axis(format="e", title="")
        ),
        alt.Color("Quantity"),
    )
)
```

```python
relax_chart.transform_filter(
    alt.FieldOneOfPredicate(field="Quantity", oneOf=["tau_D1_mean", "inv_diffusion"])
)
```

```python
long_relax.Quantity.unique()
```

```python
relax_chart.transform_filter(
    alt.FieldOneOfPredicate(
        field="Quantity", oneOf=["tau_D", "tau_Fmol", "tau_S", "tau_F"]
    )
)
```

```python
relax_chart.transform_filter(
    alt.FieldOneOfPredicate(field="Quantity", oneOf=["tau_R1", "tau_T2"])
)
```

```python
relax_chart.transform_filter(
    alt.FieldOneOfPredicate(field="Quantity", oneOf=["tau_R2", "tau_T4"])
)
```

```python
relax_chart.transform_filter(
    alt.FieldOneOfPredicate(field="Quantity", oneOf=["tau_R2", "inv_diffusion"])
)
```

```python
relax_chart.transform_filter(
    alt.FieldOneOfPredicate(field="Quantity", oneOf=["tau_D", "tau_Fmol"])
)
```
