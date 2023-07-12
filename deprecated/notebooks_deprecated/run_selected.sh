#!/bin/bash

# Make sure that conda env is activated, and screen session is started.


# human=(human_*evaluation*)
# yeast=(yeast_*evaluation*)
# files=("${human[@]}" "${human[@]}")
files=( 1_swissprot_update.ipynb 2_swissprot_keywords_go_comparison.ipynb 3_trembl_analysis.ipynb )

for file in "${files[@]}"; do
    # echo "$file"
    jupyter nbconvert --to notebook --inplace --execute "$file"
done