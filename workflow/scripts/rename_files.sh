#!/bin/bash

path=$1
file_name=$2
reg=$3

sam=$(echo $reg| cut -f1 -d'.')
h=$(echo $reg| cut -f2 -d'.')
lane=$(echo $reg| cut -f3 -d'.')
sn=$(echo $reg| cut -f4 -d'.')

expr=$(echo "$lane*$h*$sam*$sn*$file_name*")
pu=$(echo $h.$lane.$sn\_)

find $path -type f -name "$expr" | while read line;
do
    name=$(echo $line | sed -E "s/(.*\/).*($sam\_).*([0-9]_*)/\1\2$pu\3/g");
    mv $line $name;
done