#!/bin/bash

path=$1
rule=$2
capture_group=$3

find $path -name "*$rule.*" | while read line; 
do
    name=$(echo $line | sed -E "s/(.*\/._?).*($capture_group)/\1\2/g");
    mv $line $name;
done
