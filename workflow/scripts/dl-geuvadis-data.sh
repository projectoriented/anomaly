#!/usr/bin/env bash

#randomly draw samples from Geuvadis RNA-seq database
#use case: ./dl-geuvadis-data.sh 15

shuf -n $1 <(cat igsr_Geuvadis.tsv | grep "HG00" | cut -f6 | sort | uniq) > 2dl.txt

cat 2dl.txt | while read line
do
    dir_path=rna-seq_samples/$line
    mkdir -p $dir_path
    # there are samples with replicates sequenced at different lab/sites, use head to just take 1 pair.
    wget --directory-prefix=$dir_path/ --input-file=<(cat igsr_Geuvadis.tsv | grep $line | cut -f1 | head -n 2)
done
