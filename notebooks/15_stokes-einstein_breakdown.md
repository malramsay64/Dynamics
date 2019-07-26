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

Stokes-Einstein Breakdown
==============

One of the major results that we observe
is what appears to be a breakdown
in the Stokes-Einstein relationship,
most likely brought about by the reversibility of
the motions that take place.
As a result of this it would be nice to have
a figure that demonstrates this breakdown
with reference to individual particle motions,
since it is the motion of individual particles
that are the cause of this breakdown.

As an individual particle property,
we really only want to follow a smallish number of particles,
probably from a single trajectory.
To have the largest likelihood of
finding something interesting,
I will be using the dynamics from T=1.30, P=13.50.


Data Ingest/Calculation
----------------------

I don't have the translations and rotations from
every timepoint saved anywhere,
so I will have to generate all the data
from the trajectories.

```python
from pathlib import Path

import gsd.hoomd
import pandas
import numpy as np

import matplotlib.pyplot as plt
import altair as alt

from sdanalysis.dynamics import Dynamics, MolecularRelaxation, LastMolecularRelaxation
from sdanalysis.StepSize import GenerateStepSeries
import sdanalysis

import sys
sys.path.append("../src")
import figures

figures.use_my_theme()

alt.data_transformers.enable("csv")
```

### Loading Data

There are a few differences with the data used for this analysis
due to the different nature of the figures presented.

The file that we are interested in is the dump file which has linear timesteps
rather than the exponential series in the trajectory files.
The reason for the linear series is that we are interested in
following the motion of the particle through many environments,
requiring short steps for all the different environments.
Additionally there is no complication with working out
which of the series that each configuration belongs to.

Along with the trajectory,
I am loading up the dynamics computed previously.
However rather than all the dynamics,
I am only loading up those from the first keyframe.
This is because these are the dynamic quantities that
correspond to the data loaded from the trajectory.

```python
pressure = "13.50"
temperature = "1.35"
directory = Path("../data/simulations/trimer/output/")
# Load up a linear trajectory for ease of analysis
file = directory / f"trajectory-Trimer-P{pressure}-T{temperature}.gsd"
stem = f"Trimer-13.50-{temperature}"
# # Read in pre-computed dynamics properties for a single trajectory
# dynamics_df = pandas.read_hdf(directory / 'dynamics' / f'trajectory-13.50-{temperature}.hdf5',
#                               'relaxations').loc[0, :]
# # Find the fastest particles
# fastest = dynamics_df.sort_values('tau_D1').head(5).index.values
```

### Computation of Dynamics

The values of the dynamics we have calculated previously have been using an exponential timescale, which while great for the more traditional dynamics quantities becomes unhelpful in computing the molecular quantities and the analysis that will be presented in the rest of this document.

I have stored the configuration at a series of linear steps and so the below analysis computes the molecular relaxations at a single temperature (1.30) using this linear scale, while also computing the displacement and rotation at each point in time allowing for the querying of the resulting dataset.

The data is stored in a sqlite database as this is a highly queryable format while still remaining compatible with a range of systems.

```python
max_steps = None
recompute = False
# Uncomment to recompute quantities this takes a long time to compute
recompute = True
if recompute:
    dataframes = []

    with gsd.hoomd.open(str(file)) as trj:
        # Read the first frame of the snapshot to set up dynamics
        snap = sdanalysis.HoomdFrame(trj[0])
        num_mols = snap.num_mols
        dyn = Dynamics(
            timestep=snap.timestep,
            box=snap.box,
            position=snap.position,
            orientation=snap.orientation,
        )
        first_passage = MolecularRelaxation(num_mols, 0.4)
        last_passage = LastMolecularRelaxation(num_mols, 0.4, 1.0)
        exists = False
        for snap in trj:
            snap = sdanalysis.HoomdFrame(snap)
            delta_time = dyn.compute_time_delta(snap.timestep)
            delta_disp = dyn.get_displacements(snap.position)
            dataframes.append(
                {
                    "mol_index": np.arange(snap.num_mols),
                    "time": delta_time,
                    "displacement": delta_disp,
                    "rotation": dyn.get_rotations(snap.orientation),
                }
            )
            exists = True
            first_passage.add(delta_time, delta_disp)
            last_passage.add(delta_time, delta_disp)
            if max_steps and snap.configuration.step > max_steps:
                break
    pdf = pandas.DataFrame(
        {
            "first_passage": first_passage.get_status(),
            "last_passage": last_passage.get_status(),
        }
    )
    # Write data to files
    df = pandas.concat([pandas.DataFrame(data) for data in dataframes])
    df.to_hdf(f"{stem}.h5", "single_motion")
    pdf.to_hdf(f"{stem}.h5", "passage_times")
else:
    df = pandas.read_hdf(f"{stem}.h5", "single_motion")
    pdf = pandas.read_hdf(f"{stem}.h5", "passage_times")
```

```python
df.set_index(["mol_index", "time"], inplace=True, drop=False)
df.info()
```

```python
# Molecules for selection
selected = [120, 130, 140]
df_fast = df.loc[selected]
```

```python
df_fast.head()
```

```python
df_fast = df_fast.set_index(pandas.TimedeltaIndex(df_fast["time"], unit="ns"))
df_fast_groups = df_fast.groupby("mol_index").resample("1ms")
```

```python
alt.Chart(df_fast_groups.mean().reset_index(drop=True)).mark_line().encode(
    x="time", y="displacement", color="mol_index:N"
)
```

```python
mol_id = 140

df_path = df_fast.loc[df_fast.mol_index.values == mol_id, :].copy()
df_path = df_path[1:]
df_path.loc[:, "time"] = np.log(df_path.loc[:, "time"].values)
df_path = df_path[:80000]
```

```python
alt.Chart(df_path).mark_point(opacity=0.4, filled=True).encode(
    x="displacement", y="rotation", color="time"
)
```

```python
df_path.columns
```

```python
alt.Chart.mar
```

```python
df_path
```

```python
pdf.sort_values("last_passage")
pdf = pdf[pdf.last_passage != 2 ** 32 - 1]
```

```python
fig, ax = plt.subplots()
ax.plot(pdf.sort_values("first_passage").reset_index(drop=True).first_passage, ".")
ax.plot(pdf.sort_values("first_passage").reset_index(drop=True).last_passage, ".")
ax.legend()
ax.set_yscale("linear")
fig.savefig("passage_times.pdf")
```

```python
pdf.sort_values("last_passage").reset_index(drop=True).plot(logy=True)
```

```python
df.info()
```

```python
rotations_first = []
rotations_last = []
for index, first, last in pdf.itertuples():
    if last == 2 ** 32 - 1:
        continue
    rotations_first.append(df.loc[index, first].rotation)
    rotations_last.append(df.loc[index, last].rotation)
```

```python
rot_first = np.array(rotations_first)
rot_last = np.array(rotations_last)
rotations = (
    pandas.DataFrame({"first passage": rot_first})
    .melt()
    .append(pandas.DataFrame({"last passage": rot_last}).melt())
)
```

```python
altair.Chart(rotations).mark_bar(opacity=0.8).encode(
    x=altair.X("value:Q", title="Angular Displacement", bin=altair.Bin(maxbins=20)),
    y=altair.Y("count():Q", title="Occurence"),
    color=altair.Color("variable", title="Time"),
)
```

```python
fig, ax = plt.subplots()
ax.hist(rot_last, range=(0, np.pi), density=True, alpha=0.7, label="Last Passage")
ax.hist(rot_first, range=(0, np.pi), density=True, alpha=0.7, label="First Passage")
ax.legend(frameon=False)
ax.set_xlabel(r"$\Delta\theta$")
fig.savefig(f"Rotational_disp-{stem}.pdf")
```

```python
```
