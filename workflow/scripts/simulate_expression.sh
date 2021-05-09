#!/usr/bin/env bash

# I might implement getopts in the future but for now, here's a workaround.
# Please make sure the genes are in caps, and you are in a working directory.
# An example of arguments to pass: ./simulate_expression.sh *.gtf *.bam TAZ MED17 PRPF31

if [ "$#" -ne 5 ]
then
  echo "Please insert 5 arguments: 1 gtf, 1 bam, 3 genes and all delimited by space"
  exit
else
  echo -e "\c"
fi

gtf=""
bam=""
genes=()
for argument in "$@"
do
  if [[ $argument == *".bam" ]]
  then
    bam=$argument
  elif [[ $argument == *".gtf" ]]
  then
    gtf=$argument
  else
    genes+=($argument)
  fi
done

sample=$(echo "${bam##*/}" | sed "s/_.*//")

for gene in "${genes[@]}"
do
  grep "$gene" "$gtf" | grep -P "\sgene\s" > "query_$gene.gtf"
  bedtools intersect -a "$bam" -b "query_$gene.gtf" > "$sample"_sliced-"$gene".bam
  samtools view -h -s .20 -@ 4 "$sample"_sliced-"$gene".bam -o "$sample"_sliced-"$gene"_d20.bam # subsample, d20 = down sampled to 20%
done

echo "now splicing out query_"*.gtf "from $bam"
bedtools intersect -v -a "$bam" -b query_*.gtf > "$sample"_dropped3.bam # splice out the 3 genes from original bam
echo "now merging them"
samtools merge -c -p -@ 4 "$sample"_concat.bam "$sample"_sliced-*_d20.bam "$sample"_dropped2.bam