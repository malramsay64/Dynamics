#
# Makefile
# Malcolm Ramsay, 2018-03-21 13:06
#

simulation_dir = data/simulations/trimer/output
analysis_dir = data/analysis/dynamics

trajectories = $(wildcard $(simulation_dir)/trajectory-Trimer-*.gsd)
analysis = $(addprefix $(analysis_dir)/, $(notdir $(trajectories:.gsd=.h5)))

#
# Running Simulations
#

.PHONY: experiment
experiment: ## Run the trimer experiment making use of already generated files where possible
	experi --use-dependencies --input-file data/simulations/trimer/experiment.yml

#
# Dynamics Analysis
#
dynamics = data/analysis/dynamics.h5
dynamics_clean = data/analysis/dynamics_clean.h5
dynamics_agg = data/analysis/dynamics_clean_agg.h5

dynamics: | $(dynamics_agg) ## Compute dynamics quantities for all parameters of the trimer molecule

$(dynamics_agg): $(dynamics_clean)
	python src/calc_dynamics.py bootstrap $<

$(dynamics_clean): $(dynamics)
	python src/calc_dynamics.py clean --min-samples 50 $<

$(dynamics): $(analysis)
	python3 src/calc_dynamics.py collate $@ $^

$(analysis_dir)/trajectory-Trimer-P1.00-%.h5: $(simulation_dir)/trajectory-Trimer-P1.00-%.gsd
	sdanalysis --keyframe-interval 200_000 --linear-steps 100 --wave-number 2.80 comp-dynamics $< $@

$(analysis_dir)/trajectory-Trimer-P13.50-%.h5: $(simulation_dir)/trajectory-Trimer-P13.50-%.gsd
	sdanalysis --keyframe-interval 200_000 --linear-steps 100 --wave-number 2.90 comp-dynamics $< $@

#
# Figures and Notebooks
#

.PHONY: figures
figures:  | ## Create all publication figures.
	# jupyter nbconvert --ExecutePreprocessor.timeout=600 --execute notebooks/20_Publication_Figures.ipynb
	# jupyter nbconvert --ExecutePreprocessor.timeout=600 --execute notebooks/21_phd_dynamics.ipynb
	python3 src/figures.py plot-rdf --num-frames 100 data/simulations/trimer/output/dump-Trimer-P13.50-T1.50.gsd figures/thesis/radial_distribution.pdf
	python3 src/figures.py plot-ssf --num-frames 100 data/simulations/trimer/output/dump-Trimer-P13.50-T1.50.gsd figures/thesis/static_structure_factor.pdf

all_notebooks = $(wildcard notebooks/*.md)

notebooks: $(all_notebooks:.md=.ipynb)

notebooks/%.ipynb: notebooks/%.md
	cd notebooks && jupytext --to ipynb --execute $(notdir $<)

clean: ## Remove generated dynamics files and revert figures to latest committed version
	rm data/analysis/dynamics_clean.h5
	git checkout figures/*

.DEFAULT_GOAL := help
.PHONY: help
help:
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF}' $(MAKEFILE_LIST)

# vim:ft=make
#
