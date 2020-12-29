#!/usr/bin/env bash

path=$1
rule=$2
lane=$3
capture_group=$4

find $path -type f -name "$lane\_*$capture_group\_*$rule*" | while read line; 
do
    name=$(echo $line | sed -E "s/(.*\/$lane.).*($capture_group)/\1\2/g");
    mv $line $name;
done
