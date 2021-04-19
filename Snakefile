import os
import json

with open("config.json", "r") as f:
	config = json.load(f)
resultsdir = config['resultsdir']
os.makedirs(resultsdir, exist_ok=True)


rule all:
	input: 
		"docs/analysis_sim_model1.html",
		expand("{resultsdir}/result_model1.rdata", resultsdir=resultsdir)

rule test_dynamics_model1:
	input:
		"scripts/dynamics_model1.r",
		"docs/test_dynamics_model1.rmd"
	output:
		"docs/test_dynamics_model1.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_dynamics_model1.rmd\", output_format=\"all\")'"

rule test_starting_conditions_model1:
	input:
		"docs/test_starting_conditions_model1.rmd"
	output:
		"docs/test_starting_conditions_model1.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_starting_conditions_model1.rmd\", output_format=\"all\")'"

rule sim_model1:
	input:
		"docs/test_dynamics_model1.html",
		"docs/test_starting_conditions_model1.html",
		"scripts/sim_model1.r",
		"scripts/dynamics_model1.r"
	output:
		expand("{resultsdir}/result_model1.rdata", resultsdir=resultsdir)
	shell:
		"cd scripts; Rscript sim_model1.r {output}"

rule analysis_sim_model1:
	input:
		expand("{resultsdir}/result_model1.rdata", resultsdir=resultsdir),
		"docs/analysis_sim_model1.rmd"
	output:
		"docs/analysis_sim_model1.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"analysis_sim_model1.rmd\", output_format=\"all\")'"


