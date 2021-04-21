#!/usr/bin/env bash

# helper script that parses slurm output for the job ID,
# and feeds it to back to snakemake/slurm for dependencies.
# This is required when you want to use the snakemake --immediate-submit option

if [[ "Submitted batch job" =~ "$@" ]]; then
  echo -n ""
else
  deplist=$(grep -Eo "\<[0-9]{7}\>" <<< "$@" | uniq | tr '\n' ',' | sed 's/.$//')
  echo -n "--dependency=afterany:$deplist"
fi;
