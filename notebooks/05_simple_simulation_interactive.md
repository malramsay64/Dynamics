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

# Simple simulation

This is a notebook containing a simple simulation
allowing this to be used to assist in the debugging of issues,
or alternatively exploring the abilities of hoomd.

```python
import hoomd
import hoomd.md
from sdrun import molecules
```

This sets up an initial condition which is very similar to that produced using sdrun,
however this is also simple to see what all the components are and modify them.

```python
context = hoomd.context.initialize("")
trimer = molecules.Trimer()
cell_len = trimer.compute_size()
uc = hoomd.lattice.unitcell(
    N=1,
    a1=[cell_len, 0, 0],
    a2=[0, cell_len, 0],
    a3=[0, 0, 1],
    dimensions=2,
    position=[[0, 0, 0]],
    type_name=["R"],
    mass=[1.0],
    moment_inertia=[trimer.moment_inertia],
)
trimer.moment_inertia
system = hoomd.init.create_lattice(uc, 5)
for particle_type in trimer.get_types():
    if particle_type != "R":
        system.particles.types.add(particle_type)

trimer.define_potential()
rigid = trimer.define_rigid()
rigid.create_bodies()

hoomd.md.integrate.mode_standard(dt=0.005)
T = hoomd.variant.linear_interp([(0, 0.1), (100, 1), (200, 1)], zero="now")

hoomd.md.integrate.npt(hoomd.group.rigid_center(), kT=T, tau=1, P=1, tauP=1)
```

To run the simulation, the `hoomd.run` function needs to be used.

```python
hoomd.comm.get_num_ranks()
```

```python
hoomd.run(100)
```

The state of the system can be saved using `system.take_snapshot()`
which then allows the evolution of the system to be tracked and compared.

```python
snapshot = system.take_snapshot()
```
