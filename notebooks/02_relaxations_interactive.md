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

# Relaxations Summary

This is an interactive summary of the relaxation quantities calculated,
intended as a means to check that the calculated quantities are sensible
before further analysis is performed.

```python
import pandas as pd
import altair as alt
import ipywidgets as widgets

import sys
sys.path.append("../src")
import figures

figures.use_my_theme()

alt.data_transformers.enable("csv")
```

<!-- #region -->
This reads the relaxation data from the dynamics.h5 file.
The relaxation data is generated by running the command

```sh
make relaxation
```

in the root of the project directory.
The calculation of the relaxations is dependent on
the calculation of the dynamics.
<!-- #endregion -->

```python
df = pd.read_hdf("../data/analysis/dynamics_clean_agg.h5", "relaxations")
```

This is the interactive figure which provides allows selecting the relaxation quantity of interest.

```python
metadata_cols = ["temperature", "pressure"]
props = widgets.ToggleButtons(
    description="Relaxation",
    options=list(set([col.split("_")[0] for col in df.columns if col not in metadata_cols])),
)


@widgets.interact(prop=props)
def create_chart(prop):
    return figures.plot_relaxations(df, prop)
```

```python
df
```

```python
df = pd.read_hdf(
    "../data/analysis/dynamics/trajectory-Trimer-P1.00-T0.50.h5",
    "molecular_relaxations",
)
```

```python
df.reset_index().groupby("level_0").mean()
```

```python
df.index.names = ("keyframe", "molecule")
```

```python
df.head()
```