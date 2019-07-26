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

```python
from pathlib import Path

import numpy as np
import pandas
import scipy.stats
import matplotlib.pyplot as plt
from cycler import cycler
```

```python
save_figures = False
```

```python
# Configure plot look

colour_cycle = plt.rcParamsDefault["axes.prop_cycle"]
line_cycle = cycler("linestyle", [":"] + ["-"] * 4 + ["--"] + ["-"] * 4)
marker_cycle = cycler("marker", ["o", "^", "D", "X", "P"] * 2)
plt.rcParams["axes.prop_cycle"] = colour_cycle + line_cycle + marker_cycle

plt.rcParams.update({"font.size": 14})
# plt.rcParams['figure.figsize'] = [600, 400]
```

```python
source_file = Path("../data/analysis/dynamics_clean.h5")

# Read in the dynamics
dynamics = pandas.read_hdf(source_file, "dynamics")
dynamics = dynamics.query("pressure == 13.50")
dynamics = dynamics.groupby(["temperature", "time"]).mean().reset_index()
dynamics.temperature = dynamics.temperature.astype(float)
```

```python
# Read in the traditional relaxation constants computed in relaxations.
relaxations = pandas.read_hdf(source_file, "relaxations")
relaxations = relaxations.query("pressure == 13.50")

relaxations.reset_index(inplace=True)
relaxations.temperature = relaxations.temperature.astype(float)
relaxations["inv_temp"] = 1 / relaxations.temperature

# Drop lowest temperature points -> not equilibrated
relaxations = relaxations.query("temperature >= 1.30 and temperature < 3.00")
dynamics = dynamics.query("temperature >= 1.30 and temperature < 3.00")
```

```python
relaxations.columns
```

## Figure 2

Changes in liquid dynamics as measured by the structural relaxation time.

The structural relaxation time is measured in two separate ways. 

- The quantity $\tau_s$ is the value for which the function $F(\tau_s) = 1/e$, where $F(t)$ is the fraction of particles whose center of mass remain within a distance of 0.4 of their initial positions
- The quantitiy $\tau_f$ which is the time taken for an individual particle to move a distance of 0.4 from it's initial position. This quantity is typically expressed as either the mean $\langle \tau_f \rangle$ or the harmonic mean $\langle 1/\tau_f \rangle$

```python
# Figure 2a

fig, ax = plt.subplots()
for label, group in dynamics.groupby("temperature"):
    if int(label * 100) % 5 == 0:
        ax.plot("time", "com_struct", marker="", label=label, data=group)
ax.set_xscale("log")
ax.set_xlabel("t")
ax.set_ylabel("F(t)")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.legend(loc="right", bbox_to_anchor=(1.3, 0.5), frameon=False)
ax.axhline(1 / np.e, linewidth=1)
if save_figures:
    fig.savefig("../figures/figure2a.pdf")
```

```python
# Figure 2b

fig, ax = plt.subplots()
ax.plot("inv_temp", "tau_S", linestyle="", data=relaxations, label=r"$\tau_s$")
ax.plot("inv_temp", "tau_D04_mean", linestyle="", data=relaxations, label=r"$\tau_f$")
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_D04_mean / relaxations.tau_D04_hmean,
    linestyle="",
    label=r"$\langle \tau_f \rangle \langle 1/\tau_f \rangle$",
)
ax.set_yscale("log")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.legend(frameon=False)
ax.set_xlabel("1/T")
if save_figures:
    fig.savefig("../figures/figure2b.pdf")
```

```python
# Figure 3

fig, ax = plt.subplots()
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_DL04_mean / relaxations.tau_D04_mean,
    label=r"$\langle \tau_L \rangle / \langle \tau_f \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_DL04_mean / relaxations.tau_DL04_hmean,
    label=r"$\langle \tau_L \rangle \langle 1/\tau_L \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_S / relaxations.tau_DL04_mean,
    label=r"$\langle \tau_s \rangle / \langle \tau_L \rangle$",
    linestyle="",
)
ax.set_yscale("log")
ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.6))
ax.set_xlabel("1/T")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure3.pdf")
```

```python
# Figure 4a

fig, ax = plt.subplots()
for label, group in dynamics.groupby("temperature"):
    if int(label * 100) % 5 == 0:
        ax.plot("time", "rot1", marker="", label=label, data=group)
ax.set_xscale("log")
ax.set_xlabel("t")
ax.set_ylabel(r"$R_1(t)$")
ax.legend(loc="right", bbox_to_anchor=(1.3, 0.5), frameon=False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.axhline(1 / np.e, lw=1)
if save_figures:
    fig.savefig("../figures/figure4a.pdf")
```

```python
# Figure 4b

fig, ax = plt.subplots()
for label, group in dynamics.groupby("temperature"):
    if int(label * 100) % 5 == 0:
        ax.plot("time", "rot2", marker="", label=label, data=group)
ax.set_xscale("log")
ax.set_xlabel("t")
ax.set_ylabel(r"$R_2(t)$")
ax.legend(loc="right", bbox_to_anchor=(1.3, 0.5), frameon=False)
ax.axhline(1 / np.e, lw=1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure4b.pdf")
```

```python
# Figure 5a

fig, ax = plt.subplots()
ax.plot("inv_temp", "tau_R1", linestyle="", data=relaxations, label=r"$\tau_{R1}$")
ax.plot("inv_temp", "tau_R2", linestyle="", data=relaxations, label=r"$\tau_{R2}$")
ax.plot("inv_temp", "tau_S", linestyle="", data=relaxations, label=r"$\tau_S$")
ax.set_yscale("log")
ax.legend(frameon=False)
ax.set_xlabel("1/T")

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure5a.pdf")
```

```python
# Figure 5b

fig, ax = plt.subplots()
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R1 / relaxations.tau_T2_mean,
    label=r"$\tau_1 / \langle \tau_{1,f} \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_T2_mean / relaxations.tau_T2_hmean,
    label=r"$\langle \tau_{1,f} \rangle \langle 1/\tau_{1,f} \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R1 / relaxations.tau_R2,
    label=r"$\tau_1 / \tau_2 $",
    linestyle="",
)
ax.set_yscale("log")
ax.legend(frameon=False)
ax.set_xlabel("1/T")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure5b.pdf")
```

```python
# Figure 5c

fig, ax = plt.subplots()
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R2 / relaxations.tau_T4_mean,
    label=r"$\tau_2 / \langle \tau_{2,f} \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_T4_mean / relaxations.tau_T4_hmean,
    label=r"$\langle \tau_{2,f} \rangle \langle 1/\tau_{2,f} \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R1 / relaxations.tau_R2,
    label=r"$\tau_1 / \tau_2 $",
    linestyle="",
)
ax.set_yscale("log")
ax.legend(frameon=False)
ax.set_xlabel("1/T")
from matplotlib.ticker import FormatStrFormatter

ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.yaxis.set_minor_formatter(FormatStrFormatter("%g"))
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure5c.pdf")
```

```python
# Fig 6

fig, ax = plt.subplots()
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R1 / relaxations.tau_S,
    label=r"$\tau_1 / \tau_{S} $",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R2 / relaxations.tau_S,
    label=r"$\tau_2 / \tau_S$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R1 / relaxations.tau_DL04_mean,
    label=r"$ \tau_1 / \langle \tau_L \rangle$",
    linestyle="",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.tau_R2 / relaxations.tau_DL04_mean,
    label=r"$ \tau_2 / \langle \tau_L \rangle$",
    linestyle="",
)
ax.legend(frameon=False)
ax.set_xlabel("1/T")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure6a.pdf")
```

```python
# Figure 6b

select_columns = ["temperature", "tau_DL04", "tau_T2", "tau_T4"]
df = pandas.read_hdf(source_file, "molecular_relaxations")[select_columns]
df.replace(2 ** 32 - 1, np.nan, inplace=True)
df.dropna(inplace=True)

T2_on_last = np.log10(df.tau_T2.values / df.tau_DL04.values)
T4_on_last = np.log10(df.tau_T4.values / df.tau_DL04.values)

del df
```

```python
fig, ax = plt.subplots()
ax.hist(T2_on_last, bins=50, density=True)
ax.set_xlabel(r"$\log(\tau_{1,f}/\tau_L)$")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xlim(-5, 5)
if save_figures:
    fig.savefig("../figures/figure6b.pdf")
```

```python
fig, ax = plt.subplots()
plt.hist(T4_on_last, 50, density=True)
ax.set_xlabel(r"$\log(\tau_{2,f}/\tau_L)$")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xlim(-5, 5)
if save_figures:
    fig.savefig("../figures/figure6c.pdf")
```

```python
# Figure 7

fig, ax = plt.subplots()
for label, group in dynamics.groupby("temperature"):
    if int(label * 100) % 5 == 0:
        ax.plot("time", "msd", marker="", label=label, data=group)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel("t")
ax.set_ylabel(r"$\langle \Delta r^2 \rangle$")
ax.legend(loc="right", bbox_to_anchor=(1.3, 0.5), frameon=False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure7.pdf")
```

```python
# Figure 8a

fig, ax = plt.subplots()
ax.plot("inv_temp", "diffusion_constant", linestyle="", data=relaxations, label=r"$D$")
ax.plot(relaxations.inv_temp, 1 / relaxations.tau_S, linestyle="", label=r"$1/\tau_S$")
ax.set_yscale("log")
ax.legend(frameon=False)
ax.set_xlabel("1/T")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure8a.pdf")
```

```python
# Figure 8b

fig, ax = plt.subplots()
ax.plot(
    relaxations.inv_temp,
    relaxations.diffusion_constant * relaxations.tau_S,
    linestyle="",
    label=r"$D\tau_S$",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.diffusion_constant * relaxations.tau_R1,
    linestyle="",
    label=r"$D\tau_1$",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.diffusion_constant * relaxations.tau_R2,
    linestyle="",
    label=r"$D\tau_2$",
)
# ax.set_yscale('log')
ax.legend(frameon=False)
ax.set_xlabel("1/T")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure8b.pdf")
```

```python
# Figure 9

tau_DL04 = pandas.read_hdf(source_file, "molecular_relaxations")
columns = ["temperature", "pressure", "tau_DL04"]
tau_DL04 = tau_DL04[columns]
tau_DL04 = tau_DL04.query("pressure == 13.50")
tau_DL04.replace(2 ** 32 - 1, np.nan, inplace=True)
tau_DL04.dropna(inplace=True)
tau_DL04.temperature = tau_DL04.temperature.astype(float)
tau_DL04 = tau_DL04.query("temperature >= 1.30 and temperature < 3.00")

values = []
random_repetitions = 10_000


def select_from_values(values, num_selected, repetitions):
    selected_values = []
    for _ in range(repetitions):
        selected_values.append(np.random.choice(values, size=num_selected).mean())
    return scipy.stats.hmean(np.array(selected_values))


for label, group in tau_DL04.groupby(["temperature", "pressure"]):
    temperature, pressure = label
    for i in [0, 1, 2, 5, 10]:
        if i == 0:
            values.append(
                {
                    "inv_temp": 1 / temperature,
                    "num_picked": i,
                    "value": scipy.stats.hmean(group.tau_DL04.values),
                }
            )
        else:
            values.append(
                {
                    "inv_temp": 1 / temperature,
                    "num_picked": i,
                    "value": select_from_values(
                        group.tau_DL04.values, i, random_repetitions
                    ),
                }
            )

# del tau_DL04
df = pandas.DataFrame.from_records(values)
```

```python
tau_DL04.isna().sum()
```

```python
tau_DL04.groupby(["temperature", "pressure"]).filter(
    lambda x: x.isna().sum().sum() == 0
)
```

```python
test_df = tau_DL04.groupby(["temperature", "pressure"]).agg(
    [scipy.stats.gmean, scipy.stats.hmean]
)
```

```python
test_df.tau_DL04.gmean / test_df.tau_DL04.hmean
```

```python
test_df
```

```python
relaxations.tau_DL04_hmean
```

```python
fig, ax = plt.subplots()
fig.frameon = False
for label, group in df.groupby("num_picked"):
    if label == 0:
        plot_label = r"$\langle \tau_L \rangle \langle 1 / \tau_L \rangle$"
    else:
        plot_label = rf"$\langle \tau_L \rangle \langle 1 / \langle \tau_L \rangle_{{{label}}} \rangle$"
    ax.plot(
        relaxations.inv_temp,
        relaxations.tau_DL04_mean.values / group.value.values,
        linestyle="",
        label=plot_label,
    )
ax.set_yscale("log")
ax.legend(frameon=False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure9.pdf")
```

```python
# Figure 10

fig, ax = plt.subplots()
fig.frameon = False
for label, group in dynamics.groupby("temperature"):
    if int(label * 100) % 5 == 0:
        ax.plot("time", "gamma", marker="", label=label, data=group)
ax.set_xscale("log")
ax.set_xlabel("t")
ax.set_ylabel(r"$\gamma$")
ax.legend(loc="right", bbox_to_anchor=(1.21, 0.5), frameon=False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure10.pdf")
```

```python
# Figure 11

fig, ax = plt.subplots()
ax.plot(
    relaxations.inv_temp,
    relaxations.max_gamma_time / relaxations.tau_D04_mean,
    linestyle="",
    label=r"$\tau_\gamma / \langle \tau_f \rangle$",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.max_gamma_time / relaxations.tau_DL04_mean,
    linestyle="",
    label=r"$\tau_\gamma/\langle\tau_L\rangle$",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.max_gamma_time / relaxations.tau_R1,
    linestyle="",
    label=r"$\tau_\gamma/\tau_1$",
)
ax.plot(
    relaxations.inv_temp,
    relaxations.max_gamma_time / relaxations.tau_R2,
    linestyle="",
    label=r"$\tau_\gamma/\tau_2$",
)
# ax.set_yscale('log')
ax.legend(frameon=False)
ax.set_xlabel("1/T")

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
if save_figures:
    fig.savefig("../figures/figure11.pdf")
```

```python
```
