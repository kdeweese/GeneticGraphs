#!/bin/bash
n=$2
m=$3
startconnected=$4
minweight=$5
maxweight=$6
count=20000
idx=0
path=/home/kdeweese/geneticbase

while [  $idx -lt $count ];
do
    ID=$1/$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1).mtx
    while [ -f $ID ];
    do
        ID=$1/$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1).mtx
    done

    $path/evolution/create_random $n $m $ID $startconnected $minweight $maxweight
    let idx=idx+1
done
