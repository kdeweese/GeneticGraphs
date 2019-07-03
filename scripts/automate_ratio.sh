#!/bin/bash
path=/home/kdeweese/geneticbase
targetpath=$1
iters=$2
method1=$3
method2=$4
mkdir $targetpath/survived
roundstart=$6
let "roundend = $roundstart+199"

if [ $7 -eq 1 ]
then
    echo "initial round"
    $path/belos/belos_kosz_evolve_ratio.exe $targetpath/testfiles $targetpath/survived $method1 $method2 $targetpath/history.txt 0 --nullfile=ones --tol=1e-12 --max-iters=$iters --belos-xml=$path/solvers/XML/belosdefault.xml --jacobi-xml=$path/solvers/XML/jacobidefault.xml --ilut-xml=$path/solvers/XML/ilutdefault.xml --rhsfile=$path/solvers/RHS/$5.mtx
    rm $targetpath/testfiles/*
    mv $targetpath/survived/* $targetpath/testfiles/
fi
for ((i=$roundstart;i<=$roundend;i++))
do
    echo "round $i"
    $path/belos/belos_kosz_evolve_ratio.exe $targetpath/testfiles $targetpath/survived $method1 $method2 $targetpath/history.txt 5 --nullfile=ones --tol=1e-12 --max-iters=$iters --belos-xml=$path/solvers/XML/belosdefault.xml --jacobi-xml=$path/solvers/XML/jacobidefault.xml --ilut-xml=$path/solvers/XML/ilutdefault.xml --rhsfile=$path/solvers/RHS/$5.mtx
    rm $targetpath/testfiles/*
    mv $targetpath/survived/* $targetpath/testfiles/
done

