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
	dynamics_analysis bootstrap $<

$(dynamics_clean): $(dynamics)
	dynamics_analysis clean --min-samples 10 $<

$(dynamics): $(analysis)
	dynamics_analysis collate $@ $^

$(analysis_dir)/trajectory-Trimer-P1.00-%.h5: $(simulation_dir)/trajectory-Trimer-P1.00-%.gsd
	sdanalysis --keyframe-interval 1_000_000 --linear-steps 100 --wave-number 2.80 comp-dynamics --scattering-function $< $@

$(analysis_dir)/trajectory-Trimer-P13.50-%.h5: $(simulation_dir)/trajectory-Trimer-P13.50-%.gsd
	sdanalysis --keyframe-interval 1_000_000 --linear-steps 100 --wave-number 2.90 comp-dynamics --scattering-function $< $@

#
# Figures and Notebooks
#

.PHONY: figures
figures: | ## Create all publication figures.
	dynamics_figures plot-rdf --num-frames 100 data/simulations/trimer/output/dump-Trimer-P13.50-T1.50.gsd figures/radial_distribution.svg
	dynamics_figures plot-ssf --num-frames 100 data/simulations/trimer/output/dump-Trimer-P13.50-T1.50.gsd figures/static_structure_factor.svg

all_notebooks = $(wildcard notebooks/[0-9][0-9]_*.md)

notebooks: $(all_notebooks:.md=.ipynb)

.PHONY: sync
sync: ## Synchronise notebook and markdown representations
	jupytext --set-formats ipynb,md $(all_notebooks)
	jupytext --set-formats ipynb,md $(all_notebooks:.md=.ipynb)
	jupytext --sync --pipe-fmt py --pipe black $(all_notebooks:.md=.ipynb)

notebooks/%.ipynb: notebooks/%.md $(dynamics_agg)
	cd notebooks && jupytext --to ipynb --execute $(notdir $<)

clean: ## Remove generated dynamics files and revert figures to latest committed version
	rm data/analysis/dynamics_clean.h5
	git checkout figures/*

#
# Reports
#

# Variables
#
# These find all the different items which are going to become reports. I am converting
# all the notebooks to a pdf report, as well as markdown files in the reports directory.
#
# Additionally where there are figures which need converting from svg to pdf, this takes
# finds the targets.

report_targets = $(patsubst %.md, %.pdf, $(wildcard notebooks/*.md))
all_figures = $(wildcard figures/*.svg)

# Commands
#
# There are a couple of commands here. Generating all the figures and additionally
# creating the pdf reports.

.PHONY: reports
reports: figures notebooks $(report_targets) ## Generate pdf reports

# Conversions
#
# These are the rules to convert the different filetypes to a pdf

%.pdf: %.md
	cd $(dir $<); pandoc $(notdir $<) --pdf-engine=tectonic --filter ../src/pandoc-svg.py -o $(notdir $@)

#
# Help
#
# The magic which provides some help output when running make by itself

.DEFAULT_GOAL := help
.PHONY: help
help:
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF}' $(MAKEFILE_LIST)

# vim:ft=make
#
