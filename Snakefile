import os
import json

with open("config.json", "r") as f:
	config = json.load(f)
resultsdir = config['resultsdir']
os.makedirs(resultsdir, exist_ok=True)


rule all:
	input: 
		"docs/analysis_stst_sim_model4.html",
		expand("{resultsdir}/result_model4.rdata", resultsdir=resultsdir)

rule test_dynamics_model4:
	input:
		"scripts/dynamics_model4.r",
		"docs/test_dynamics_model4.rmd"
	output:
		"docs/test_dynamics_model4.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_dynamics_model4.rmd\", output_format=\"all\")'"

rule sim_model4:
	input:
		"docs/test_dynamics_model4.html",
		"scripts/sim_model4.r",
		"scripts/dynamics_model4.r"
	output:
		expand("{resultsdir}/result_model4.rdata", resultsdir=resultsdir)
	shell:
		"cd scripts; Rscript sim_model4.r {output}"

rule analysis_dyn_sim_model4:
	input:
		"docs/analysis_dyn_sim_model4.rmd"
	output:
		"docs/analysis_dyn_sim_model4.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"analysis_dyn_sim_model4.rmd\", output_format=\"all\")'"

rule analysis_stst_sim_model4:
	input:
		expand("{resultsdir}/result_model4.rdata", resultsdir=resultsdir),
		"docs/analysis_stst_sim_model4.rmd"
	output:
		"docs/analysis_stst_sim_model4.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"analysis_stst_sim_model4.rmd\", output_format=\"all\")'"


