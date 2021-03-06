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

# Relaxation Dynamics

There are many different methods
which can be used to understand the dynamics of a system.
From short timescale events like structural relaxation,
to long timescale events including dynamics.
There are also important degrees of freedom
in the rotations which have their own relaxation timescales.

This is a collection of figures and analysis
for the understanding of relaxation over a series of variables.

All the simulation data is from a set of simulations stored on RDS
in the folder `data/simulations/trimer`.
This set of simulations were all run at a pressure of 13.50.


## Setup

This imports the modules necessary for running the code in the rest of this notebook
while also setting up some helper functions to make the rest of the analysis simpler.

```python
import pandas as pd
import numpy as np
import altair as alt
from dynamics_analysis import figures
```


### Import Data

This reads all the data from disk,
creating a pandas DataFrame for both
the dynamics and the relaxation quantities.

```python
data_file = "../data/analysis/dynamics_clean_agg.h5"
```

```python
dynamics = pd.read_hdf(data_file, "dynamics")
dynamics = dynamics.query("pressure == 13.50")
dynamics = dynamics.query("temperature > 1.25")
dynamics = dynamics.groupby(["temperature", "pressure", "time"]).mean()
dynamics = dynamics.reset_index()
```

```python
relaxations = pd.read_hdf(data_file, "relaxations")
# relaxations = relaxations.query("pressure == 13.50")
```

```python
mol_relax = pd.read_hdf(data_file, "molecular_relaxations")
# mol_relax = mol_relax.query("pressure == 13.50")
```

The available temperatures for plotting are listed below.
To make the visualisations of time dependent properties less cluttered
only a subset of the temperatures will be used for plotting.
For relaxation timescales all the temperatures will be used.
This list of temperatures for plotting
can be modified by changing the plot_temperatures variable.

```python
np.sort(dynamics["temperature"].unique())
```

### Diffusion


Generating the figures for diffusive relaxation of molecules.
Lines are indications of fit,
where $D t$ is the function defining the fit.
The parameter `D` is used as the diffusion constant in the following figure.

```python
fig = figures.plot_dynamics(dynamics, "msd", title="Mean Square Displacement", scale="log")
x=np.linspace(1e-2, 1, 5)
df = pd.DataFrame({"x2": x, "y2": x**2, "x": x*5e6, "y": x*50})
c = alt.Chart(df).encode(x="x").mark_line(color="black")
fig + c.encode(x="x2", y="y2").mark_line(color="black") + c.encode(x="x", y="y")
```

### Rotational Relaxation R2

This is the second order relaxation function,
given by;

$$R_2(t) = \langle 2[ \hat{\mathbf{e}}(0) \cdot \hat{\mathbf{e}}(t)]^2 - 1 \rangle$$

where the value is averaged over all molecules and starting configurations.

```python
figures.plot_dynamics(dynamics, "rot2", title="Rotational Relaxation")
```

### Structural Relaxation


This is the fraction of particles which
have moved a distance of 0.3 from their initial position.

```python
figures.plot_dynamics(dynamics, "struct", title="Intermediate Scattering Function")
```

## Summary Values


The values calculated below summarise the above information,
providing a method of investigating temperature dependence of these properties.

```python
figures.plot_relaxations(relaxations, "msd", title="Diffusion Constant")
```

```python
relaxations.columns
```

```python
figures.plot_relaxations(relaxations, "inv_diffusion", "Inverse Diffusion Constant")
```

```python
figures.plot_relaxations(relaxations, "rot2", title="Rotational Relaxation Time")
```

```python
mol_relax.columns
```

```python
figures.plot_relaxations(mol_relax, "tau_T2")
```

```python
relaxations.columns
```

```python
values = pd.DataFrame(
    {
        "inv_temp_norm": relaxations.inv_temp_norm,
        "Temperature": relaxations.temperature,
        "r1r2": relaxations.rot1_mean / relaxations.rot2_mean,
        "Dr1T": relaxations.msd_mean
        * relaxations.rot2_mean
        * relaxations.inv_temp_norm,
        "Dr2T": relaxations.msd_mean
        * relaxations.rot2_mean
        * relaxations.inv_temp_norm,
        "DsT": relaxations.msd_mean
        * relaxations.com_struct_mean
        * relaxations.inv_temp_norm,
        "Dr1": relaxations.msd_mean * relaxations.rot1_mean,
        "Dr2": relaxations.msd_mean * relaxations.rot2_mean,
    }
)
```

```python
c = (
    alt.Chart(values)
    .mark_point(filled=True)
    .encode(alt.X("inv_temp_norm", scale=alt.Scale(zero=False)))
)

c.encode(y="r1r2")
```
