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

# Dynamics Figures

These are a collection of figures for my PhD thesis.
The naming of the figures will remain consistent throughout and the ordering of the figures should be roughly the same as in the thesis.

```python
%load_ext autoreload
%autoreload 2
```

```python
# Make dealing with filesystem paths much, much simpler
from pathlib import Path

# Read/write data files and data analysis
import pandas
import numpy as np
import altair as alt

import sys
sys.path.append("../src")
import figures

figures.use_my_theme()

alt.data_transformers.enable('csv')
```

This notebook generates a collection of figures which are exported as pdf files to the `../figures/thesis` directory.
To prevent the figures from being saved, only outputting the resulting files to this notebook
you can set the `save_figures` variable to `False`.

```python
save_figures = True
save_figures = False
```

## Datasets

The datasets used for the generation of these figures are found in the `../data/analysis` directory,
with the `dynamics.h5` file containing all the dynamics results in 3 tables

- `dynamics` -> containing the standard raw dynamics quantities including the mean-squared-displacement, structural relaxation, and many others.
    These values include each of these quantities for a number of starting configurations, allowing for the calculation of errors for these quantities.
    For a full list of the available quantities see the [dynamics_interactive notebook](01_dynamics_interactive.ipynb).
- `molecular_relaxations` -> containing the molecular relaxation values for each of the molecules
    allowing comparisons between quantities for a single molecule.
- `relaxations` -> containing the aggregated relaxation value for each of the quantities in the `dynamics` and `molecular_relaxations` tables.
    All these values can be investigated in the [relaxations_interactive notebook](02_relaxations_interactive.ipynb)


```python
# Where the data files with the results are located
data_dir = Path('../data/analysis')

# Load data for most of the figures
dynamics_df = pandas.read_hdf(data_dir / 'dynamics_clean_agg.h5', 'dynamics')

dynamics_df = dynamics_df.query("pressure == 13.50")

# Output path for all figures
figure_dir = Path("../figures/thesis")
# Ensure the directory exists
figure_dir.mkdir(exist_ok=True)
```

## Normalisation by Melting Point

The dynamics of many quantities are plotted as a fraction of the melting point,
this creates a column in the data frames with this normalised temperature.

```python
dynamics_df['temp_norm'] = 0

select_high_pressure = (dynamics_df.pressure == 13.50).values
dynamics_df.loc[select_high_pressure, 'temp_norm'] = 1.35 / dynamics_df.loc[select_high_pressure, 'pressure']

select_low_pressure = (dynamics_df.pressure == 1.00).values
dynamics_df.loc[select_low_pressure, 'temp_norm'] = 0.36 / dynamics_df.loc[select_low_pressure, 'pressure']
```

```python
dynamics_df.info()
```

## Comparative Dynamics

These are a collection of dynamics quantities to establish that the system we are dealing with has behaviour that more or less aligns with much of the literature.

```python
from figures import plot_dynamics
```

### Mean Squared Displacement


```python
c = plot_dynamics(dynamics_df, 'msd', scale='log')

if save_figures:
    c.save(str(figure_dir / "mean_squared_displacement.svg"), webdriver='firefox')

c
```

## Non-gaussian

```python
c = plot_dynamics(dynamics_df, 'alpha', title="Non Gaussian")

if save_figures:
    c.save(str(figure_dir / "non_gaussian.svg"), webdriver='firefox')

c
```

## Structural Relaxation

```python
c = plot_dynamics(dynamics_df, 'struct', title="Structrual Relaxation")

if save_figures:
    c.save(str(figure_dir / "structural_relaxation.svg"), webdriver='firefox')

c
```

```python
c = plot_dynamics(dynamics_df, 'scattering_function', title="Intermediate Scattering Function")

if save_figures:
    c.save(str(figure_dir / "scattering_function.svg"), webdriver='firefox')

c
```

## Rotational Relaxation


```python
c = plot_dynamics(dynamics_df, 'rot2', title="Rotational Relaxation")

if save_figures:
    c.save(str(figure_dir / "rotational_relaxtion.svg"), webdriver='firefox')

c
```

## Relaxation Quantities

```python
relaxations_df = pandas.read_hdf(data_dir / "dynamics_clean_agg.h5", "relaxations")
# relaxations_df = relaxations_df.query("pressure == 13.50")
mask = relaxations_df.pressure == 13.50
t_melting_high = 1.35
t_melting_low = 0.35
relaxations_df['inv_temp'] = 0.
relaxations_df.loc[mask, 'inv_temp'] = t_melting_high / relaxations_df.temperature
relaxations_df.loc[~mask, 'inv_temp'] = t_melting_low / relaxations_df.temperature
relaxations_df[relaxations_df < 0] = np.NaN

# relaxations_df['inv_diffusion'] = 1 / relaxations_df.diffusion_constant
```

```python
from figures import plot_relaxations
```

```python
relaxations_df.head()
```

```python
relaxations_df.columns
```

### Scattering Function


```python
plot_relaxations(relaxations_df, 'scattering_function')
```

```python
plot_relaxations(relaxations_df, 'struct', title="Structural Relaxation")

if save_figures:
    c.save(str(figure_dir / "isf_relaxation.svg"), webdriver='firefox')

c
```

### Diffusion


```python
c = plot_relaxations(relaxations_df, 'inv_diffusion', title="1/D")

if save_figures:
    c.save(str(figure_dir / "diffusion_constant.svg"), webdriver='firefox')

c
```

### Rotational Relaxation

```python
c = plot_relaxations(relaxations_df, 'rot2', title="Rotational Relaxation")

if save_figures:
    c.save(str(figure_dir / "rotational_relaxation.svg"), webdriver='firefox')

c
```

## Molecular Relaxations

```python
mol_df = pandas.read_hdf(data_dir / "dynamics_clean_agg.h5", "molecular_relaxations")

mask = mol_df.pressure == 13.50
t_melting_high = 1.35
t_melting_low = 0.35
mol_df['inv_temp'] = 0.
mol_df.loc[mask, 'inv_temp'] = t_melting_high / mol_df.temperature
mol_df.loc[~mask, 'inv_temp'] = t_melting_low / mol_df.temperature
```

```python
relax_df = relaxations_df.set_index(['temperature', 'pressure', 'inv_temp']).join(
    mol_df.set_index(['temperature', 'pressure', 'inv_temp'])
).reset_index()
```

```python
from figures import reshape_dataframe, plot_multi_relaxations
```

```python
mol_df
```

```python
melt_df = reshape_dataframe(relax_df)
```

```python
melt_df.variable.unique()
```

```python
plot_relaxations(mol_df, "tau_F")
```

```python
(
    plot_multi_relaxations(melt_df, "tau_F", title="Molecular") +
    plot_multi_relaxations(melt_df, "scattering_function", title="Scattering")
)
```

```python
comp_relax_df.variable.unique()
```

```python
import altair as alt

c = alt.Chart(comp_relax_df).mark_point().encode(
    alt.X('inv_temp:Q', title="Tm/T"),
    alt.Y('value:Q', title="Relaxtion Time", scale=alt.Scale(type='log'), axis=alt.Axis(format='e')),
    alt.Color('pressure:N', title="Pressure"),
    alt.Shape('variable:N', title="Quantity"),
)

c.transform_filter((alt.datum.variable == "rot1_value")) + c.transform_filter((alt.datum.variable == "rot2_value"))
```

```python
import altair as alt

c = alt.Chart(comp_relax_df).mark_point().encode(
    alt.X('inv_temp', title="Tm/T"),
    alt.Y('value', title="Relaxtion Time", scale=alt.Scale(type='log'), axis=alt.Axis(format='e')),
    alt.Color('pressure:N', title="Pressure"),
    alt.Shape('variable', title="Quantity"),
)

c.transform_filter((alt.datum.variable == "inv_diffusion") | (alt.datum.variable == "tau_D"))
```

```python
import altair as alt

c = alt.Chart(comp_relax_df).mark_point().encode(
    alt.X('inv_temp', title="Tm/T"),
    alt.Y('value', title="Relaxtion Time", scale=alt.Scale(type='log'), axis=alt.Axis(format='e')),
    alt.Color('pressure:N', title="Pressure"),
    alt.Shape('variable', title="Quantity"),
)

c.transform_filter((alt.datum.variable == "inv_diffusion") | (alt.datum.variable == "tau_D"))
```

```python
c = alt.Chart(comp_relax_df).mark_point().encode(
    alt.X('inv_temp', title="Tm/T"),
    alt.Y('value', title="Relaxtion Time", scale=alt.Scale(type='log'), axis=alt.Axis(format='e')),
    alt.Color('pressure:N', title="Pressure"),
    alt.Shape('variable', title="Quantity"),
)

c.transform_filter((alt.datum.variable == "tau_Fmol"))

```

```python
fig, ax = plt.subplots()

ax.plot('inv_temp', 'inv_diffusion', 's', data=relaxations_df, label="1/D")
ax.plot('inv_temp', 'tau_R1', 'x', data=relaxations_df, label=r"$\tau_S$")
ax.plot('inv_temp', 'tau_F', 'o', data=relaxations_df, label=r"$\tau_F$")
ax.plot('inv_temp', 'tau_S', 'o', data=relaxations_df, label=r"$\tau_S$")

ax.set_xlabel('$T_m/T$')
ax.set_ylabel("Relaxation Time")

# ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

fig.savefig("Relaxations.pdf")
```

```python
c = alt.Chart(comp_relax_df).mark_point().encode(
    alt.X('inv_temp', title="Tm/T"),
    alt.Y('value', title="Relaxtion Time", scale=alt.Scale(type='log'), axis=alt.Axis(format='e')),
    alt.Color('pressure:N', title="Pressure"),
    alt.Shape('variable', title="Quantity"),
)

c.transform_filter((alt.datum.variable == "inv_diffusion"))
```

```python
dynamics_df.info()
```

```python
relaxations_df
```
