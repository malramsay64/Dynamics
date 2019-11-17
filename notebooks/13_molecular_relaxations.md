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

# Molecular Relaxations

Rather than investigating relaxations as
an average over all the molecules,
we want to investigate individual molecules
and how their behaviour is related to that of the aggregate.

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
mol_relax = pd.read_hdf(data_file, "molecular_relaxations")
relax_df = relaxations.merge(mol_relax, on=["temperature", "pressure", "inv_temp_norm"])
```

```python
plot_relax_df = figures.reshape_dataframe(relax_df)
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["tau_D", "inv_diffusion"])
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["tau_D", "inv_diffusion"])
```

```python
plot_relax_df.variable.unique()
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["tau_D", "struct", "tau_F"])
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["tau_T2", "rot1"])
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["tau_T4", "rot2"])
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["rot2", "inv_diffusion"])
```

```python
figures.plot_multi_relaxations(plot_relax_df, ["tau_D"])
```
