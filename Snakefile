import os
import json

with open("config.json", "r") as f:
	config = json.load(f)
resultsdir = config['resultsdir']
os.makedirs(resultsdir, exist_ok=True)


rule all:
	input: 
<<<<<<< HEAD
		"docs/analysis_stst_sim_model4.html",
=======
		"docs/analysis_sim_model4.html",
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
		expand("{resultsdir}/result_model4.rdata", resultsdir=resultsdir)

rule test_dynamics_model4:
	input:
		"scripts/dynamics_model4.r",
		"docs/test_dynamics_model4.rmd"
	output:
		"docs/test_dynamics_model4.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_dynamics_model4.rmd\", output_format=\"all\")'"

<<<<<<< HEAD
rule sim_model4:
	input:
		"docs/test_dynamics_model4.html",
=======
rule test_starting_conditions_model4:
	input:
		"docs/test_starting_conditions_model4.rmd"
	output:
		"docs/test_starting_conditions_model4.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"test_starting_conditions_model4.rmd\", output_format=\"all\")'"

rule sim_model4:
	input:
		"docs/test_dynamics_model4.html",
		"docs/test_starting_conditions_model4.html",
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
		"scripts/sim_model4.r",
		"scripts/dynamics_model4.r"
	output:
		expand("{resultsdir}/result_model4.rdata", resultsdir=resultsdir)
	shell:
		"cd scripts; Rscript sim_model4.r {output}"

<<<<<<< HEAD
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
=======
rule analysis_sim_model4:
	input:
		expand("{resultsdir}/result_model4.rdata", resultsdir=resultsdir),
		"docs/analysis_sim_model4.rmd"
	output:
		"docs/analysis_sim_model4.html"
	shell:
		"cd docs; Rscript -e 'rmarkdown::render(\"analysis_sim_model4.rmd\", output_format=\"all\")'"
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528


