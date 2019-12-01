---
jupyter:
  jupytext:
    formats: ipynb,md
    target_format: ipynb,md
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

# Dynamics Figures

These are a collection of figures for my PhD thesis.
The naming of the figures will remain consistent throughout and the ordering of the figures should be roughly the same as in the thesis.

```python
# Make dealing with filesystem paths much, much simpler
from pathlib import Path

# Read/write data files and data analysis
import pandas
import numpy as np
import altair as alt
from dynamics_analysis import figures, calc_dynamics, util

```

This notebook generates a collection of figures
which are exported as svg files
to the `../figures` directory.

## Datasets

The datasets used for the generation of these figures
are found in the `../data/analysis` directory,
with the `dynamics_clean_agg.h5` file containing all the dynamics results in 3 tables

- `dynamics` -> containing the standard raw dynamics quantities including the mean-squared-displacement, structural relaxation, and many others.
    These values include each of these quantities for a number of starting configurations, allowing for the calculation of errors for these quantities.
    For a full list of the available quantities see the [dynamics_interactive notebook](01_dynamics_interactive.ipynb).
- `molecular_relaxations` -> containing the molecular relaxation values for each of the molecules
    allowing comparisons between quantities for a single molecule.
- `relaxations` -> containing the aggregated relaxation value for each of the quantities in the `dynamics` and `molecular_relaxations` tables.
    All these values can be investigated in the [relaxations_interactive notebook](02_relaxations_interactive.ipynb)

```python
# Where the data files with the results are located
data_dir = Path("../data/analysis")

# Load data for most of the figures
dynamics_df = pandas.read_hdf(data_dir / "dynamics_clean_agg.h5", "dynamics")

dynamics_df = dynamics_df.query("pressure == 13.50")
dynamics_df = dynamics_df.sort_values("time")

# Output path for all figures
figure_dir = Path("../figures")
# Ensure the directory exists
figure_dir.mkdir(exist_ok=True)
```

## Comparative Dynamics

These are a collection of dynamics quantities
to establish that the system we are dealing with
has behaviour that aligns with much of the literature.

### Mean Squared Displacement

```python
c = figures.plot_dynamics(
    dynamics_df, "msd", title="Mean Squared Displacement", scale="log"
)

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "mean_squared_displacement.svg"), webdriver="firefox")
```

![mean squared displacement](../figures/mean_squared_displacement.svg)


## Non-gaussian

```python
c = figures.plot_dynamics(dynamics_df, "alpha", title="Non Gaussian")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "non_gaussian.svg"), webdriver="firefox")
```

![non-gaussian](../figures/non_gaussian.svg)


## Structural Relaxation

```python
c = figures.plot_dynamics(dynamics_df, "struct", title="Structrual Relaxation")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "structural_relaxation.svg"), webdriver="firefox")
```

![structural relaxation](../figures/structural_relaxation.svg)

```python
c = figures.plot_dynamics(
    dynamics_df, "scattering_function", title="Intermediate Scattering Function"
)

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "scattering_function.svg"), webdriver="firefox")
```

![intermediate scattering function](../figures/scattering_function.svg)

## Rotational Relaxation

```python
c = figures.plot_dynamics(dynamics_df, "rot2", title="Rotational Relaxation")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "rotational_relaxation.svg"), webdriver="firefox")
```

![rotational relaxation](../figures/rotational_relaxation.svg)

```python
c = figures.plot_dynamics(dynamics_df, "alpha_rot", title="Rotational Non-Gaussian")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "rotational_alpha.svg"), webdriver="firefox")
```

![rotational  heterogenaeties](../figures/rotational_alpha.svg)


## Rotational and Translational Coupling

```python
c = figures.plot_dynamics(
    dynamics_df, "gamma", title="Translational-Rotational Coupling"
)

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "gamma.svg"), webdriver="firefox")
```

![Tranlational and Rotational Coupling](../figures/gamma.svg)


## Relaxation Quantities

```python
relaxations_df = pandas.read_hdf(data_dir / "dynamics_clean_agg.h5", "relaxations")
# relaxations_df = relaxations_df.query("pressure == 13.50")
relaxations_df[relaxations_df < 0] = np.NaN
```

```python
relaxations_df.columns
```

### Scattering Function

```python
c = figures.plot_relaxations(relaxations_df, "scattering_function")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "scattering_function_summary.svg"), webdriver="firefox")
```

![scattering function relaxation times](../figures/scattering_function_summary.svg)

```python
c2 = figures.plot_relaxations(
    relaxations_df, "com_struct", title="Structural Relaxation"
)

with alt.data_transformers.enable("default"):
    c2.save(str(figure_dir / "structural_relaxation_summary.svg"), webdriver="firefox")
```

```python
relaxations_df
```

![structural relaxation times](../figures/structural_relaxation_summary.svg)


### Diffusion

```python
c = figures.plot_relaxations(relaxations_df, "inv_diffusion", title="1/D")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "diffusion_constant_summary.svg"), webdriver="firefox")
```

![Diffusion constant](../figures/diffusion_constant_summary.svg)


### Rotational Relaxation

```python
c = figures.plot_relaxations(relaxations_df, "rot2", title="Rotational Relaxation")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "rotational_relaxation_summary.svg"), webdriver="firefox")
```

![rotational relaxation times](../figures/rotational_relaxation_summary.svg)


## Heterogeneities


```python
c = figures.plot_relaxations(relaxations_df, "alpha", title="Rotational Relaxation")

with alt.data_transformers.enable("default"):
    c.save(str(figure_dir / "alpha_relaxation_summary.svg"), webdriver="firefox")
c
```
