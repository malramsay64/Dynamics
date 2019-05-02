#
# Makefile
# Malcolm Ramsay, 2018-03-21 13:06
#
dynamics = data/analysis/dynamics.h5
dynamics_clean = data/analysis/dynamics_clean.h5

simulation_dir = data/simulations/trimer/output
analysis_dir = data/analysis/dynamics

trajectories = $(wildcard $(simulation_dir)/trajectory-Trimer-*.gsd)
analysis = $(addprefix $(analysis_dir)/, $(notdir $(trajectories:.gsd=.h5)))


.PHONY: experiment
experiment: ## Run the trimer experiment making use of already generated files where possible
	experi --use-dependencies --input-file data/simulations/trimer/experiment.yml

${dynamics}: $(analysis)
	python3 src/collate_data.py $@ $^

$(analysis_dir)/trajectory-Trimer-P1.00-%.h5: $(simulation_dir)/trajectory-Trimer-P1.00-%.gsd
	sdanalysis --keyframe-interval 200_000 --wave-number 2.80 comp-dynamics $< $@

$(analysis_dir)/trajectory-Trimer-P13.50-%.h5: $(simulation_dir)/trajectory-Trimer-P13.50-%.gsd
	sdanalysis --keyframe-interval 200_000 --wave-number 2.90 comp-dynamics $< $@

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

clean: ## Remove generated dynamics files and revert figures to latest committed version
	rm data/analysis/dynamics_clean.h5
	git checkout figures/*

.DEFAULT_GOAL := help
.PHONY: help
help:
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF}' $(MAKEFILE_LIST)

# vim:ft=make
#
