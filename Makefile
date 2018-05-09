#
# Makefile
# Malcolm Ramsay, 2018-03-21 13:06
#

dynamics:
	ls data/simulations/trimer/trajectory-* | xargs -n1 sdanalysis comp_dynamics -o data/analysis

relaxations: data/analysis/dynamics.h5
	python src/relaxations.py

data/analysis/dynamics.h5: dynamics

.PHONY: relaxations dynamics

# vim:ft=make
#
