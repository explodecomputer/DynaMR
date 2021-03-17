import os
import json

with open("config.json", "r") as f:
	config = json.load(f)
resultsdir = config['resultsdir']
os.makedirs(resultsdir, exist_ok=True)


rule all:
	input: "docs/analysis_sim.html"


rule test_dynamics:
	input:
		"scripts/dynamics.r",
		"docs/test_dynamics.rmd"
	output:
		"docs/test_dynamics.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_dynamics.rmd\", output_format=\"all\")'"

rule test_starting_conditions:
	input:
		"docs/test_starting_conditions.rmd"
	output:
		"docs/test_starting_conditions.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_starting_conditions.rmd\", output_format=\"all\")'"

rule sim:
	input:
		"docs/test_dynamics.html",
		"docs/test_starting_conditions.html",
	output:
		expand("{resultsdir}/result.rdata", resultsdir=resultsdir)
	shell:
		"cd scripts; Rscript sim.r {output}"

rule analysis_sim:
	input:
		expand("{resultsdir}/result.rdata", resultsdir=resultsdir),
		"docs/analysis_sim.rmd"
	output:
		"docs/analysis_sim.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"analysis_sim.rmd\", output_format=\"all\")'"


