#!/bin/bash

# nohup bash fit_structure_job.sh > /dev/null 2>&1 &

for pp in {1..10}; do
    sed "s/pp = \"pathway_information\"/pp = ${pp}/" ./template/fit_structure_template.R > fit_structure_pp${pp}.R
done

for pp in {1..10}; do
    Rscript fit_structure_pp${pp}.R > ./log/output_fit_structure_pp${pp}.log 2>&1 &
    sleep 5
done

wait

for pp in {1..10}; do
    rm fit_structure_pp${pp}.R
done