#
# Makefile
# Malcolm Ramsay, 2018-03-21 13:06
#

.PHONY: dynamics
dynamics: ## Compute dynamics quantities for all parameters of the trimer molecule
	sdanalysis comp_dynamics -o data/analysis data/simulations/trimer/output/trajectory-*

.PHONY: relaxations
relaxations: ## Compute the summary relaxation timescales of the dynamic quantitites
	sdanalysis comp_relaxations data/analysis/dynamics.h5

.PHONY: figures
figures:  ## Create all publication figures.
	jupyter nbconvert --ExecutePreprocessor.timeout=600 --execute notebooks/03_Publication_Figures.ipynb

.DEFAULT_GOAL := help
.PHONY: help
help:
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF}' $(MAKEFILE_LIST)

# vim:ft=make
#
