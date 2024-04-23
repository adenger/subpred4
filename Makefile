.PHONY: env_export clear_tmp_files requirements package data_import blast_databases data_export

#################################################################################
# Setup                                                                         #
#################################################################################

## Import raw data from archive created with data_export
data_import:
	tar xvf subpred4_data.tar.gz

## Install packages with mambaforge (conda also works but takes hours)
requirements:
	# conda update -n base -c defaults conda
	mamba env create --file environment.yml

## Install code as python package to use it in notebooks
package:
	pip install -e .

## Export current env to new file
env_export:
	conda env export | sed '/^prefix/d' > environment.yml
	# conda env export --from-history | sed '/^prefix/d' > environment_history.yml

## Export raw data for sharing (requires active conda env)
data_export:
	tar -cvf - -T data_backup_list.txt | pigz > subpred4_data.tar.gz

## Create local databases for generating new PSSM files
blast_databases:
	cd data/raw/uniref/uniref50 && makeblastdb -in uniref50.fasta -parse_seqids -dbtype prot
	cd data/raw/uniref/uniref90 && makeblastdb -in uniref90.fasta -parse_seqids -dbtype prot

## Clean up tmp files that are not needed
clear_tmp_files:
	find data/intermediate/blast -name "*.log" -delete
	find data/intermediate/blast -name "*.fasta" -delete

#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Copied from https://github.com/drivendata/cookiecutter-data-science
# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
