#!/bin/bash

# nohup bash fit_outcome_job.sh > /dev/null 2>&1 &

drugs=("erlotinib" "gefitinib" "afatinib" "dacomitinib" "osimertinib" "cisplatin")

for drug in "${drugs[@]}"; do
    for pp in {1..10}; do
        sed "s/pp = \"pathway_information\"/pp = ${pp}/; s/drug_name = \"drug_name\"/drug_name = \"${drug}\"/" ./template/fit_outcome_template.R > fit_outcome_${drug}_pp${pp}.R
    done
done

for drug in "${drugs[@]}"; do
    for pp in {1..10}; do
        Rscript fit_outcome_${drug}_pp${pp}.R > ./log/output_fit_outcome_${drug}_pp${pp}.log 2>&1 &
    done
done

wait

for drug in "${drugs[@]}"; do
    for pp in {1..10}; do
        rm fit_outcome_${drug}_pp${pp}.R
    done
done