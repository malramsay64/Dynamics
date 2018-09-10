#
# Makefile
# Malcolm Ramsay, 2018-03-21 13:06
#

dynamics:
	ls data/simulations/trimer/trajectory-* | xargs -n1 sdanalysis comp_dynamics -o data/analysis

relaxations:
	sdanalysis comp_relaxations data/analysis/dynamics.h5

figures:
	(cd notebooks; jupyter nbconvert --ExecutePreprocessor.timeout=600 --execute Figures.ipynb)

.PHONY: relaxations dynamics figures

# vim:ft=make
#
