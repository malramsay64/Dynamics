#
# Makefile
# Malcolm Ramsay, 2018-03-21 13:06
#
dynamics = data/analysis/dynamics.h5
dynamics_clean = data/analysis/dynamics_clean.h5

.PHONY: experiment
experiment: ## Run the trimer experiment making use of already generated files where possible
	experi --use-dependencies --input-file data/simulations/trimer/experiment.yml

${dynamics}:
	sdanalysis --keyframe-interval 200_000 comp-dynamics -o data/analysis data/simulations/trimer/output/trajectory-*

${dynamics_clean}: ${dynamics}
	python src/data_cleanup.py $<

dynamics: | ${dynamics_clean} ## Compute dynamics quantities for all parameters of the trimer molecule

.PHONY: relaxations
relaxations: | ${dynamics_clean} ## Compute the summary relaxation timescales of the dynamic quantitites
	sdanalysis comp-relaxations ${dynamics_clean}

.PHONY: figures
figures:  | ## Create all publication figures.
	jupyter nbconvert --ExecutePreprocessor.timeout=600 --execute notebooks/20_Publication_Figures.ipynb
	jupyter nbconvert --ExecutePreprocessor.timeout=600 --execute notebooks/21_phd_dynamics.ipynb

.DEFAULT_GOAL := help
.PHONY: help
help:
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF}' $(MAKEFILE_LIST)

# vim:ft=make
#
