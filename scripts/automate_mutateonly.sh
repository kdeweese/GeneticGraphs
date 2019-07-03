#!/bin/bash
iters=$2
path=/home/kdeweese/geneticbase
targetpath=$1
method1=$3
mkdir $targetpath/survived
roundstart=$6
let "roundend = $roundstart+200"
if [ $7 -eq 1 ]
then
    echo "initial round"
    $path/solvers/belos_kosz_mutate.exe $targetpath/testfiles $targetpath/survived $method1 $4 0 --nullfile=ones --tol=1e-12 --max-iters=$iters --belos-xml=$path/solvers/XML/belosdefault.xml --jacobi-xml=$path/solvers/XML/jacobidefault.xml --ilut-xml=$path/solvers/XML/ilutdefault.xml --rhsfile=$path/solvers/RHS/$5.mtx
    rm $targetpath/testfiles/*
    #mkdir $targetpath/survivedinit/
    #cp $targetpath/survived/* $targetpath/survivedinit/
    mv $targetpath/survived/*.mtx $targetpath/testfiles/
fi
for ((i=$roundstart;i<=$roundend;i++))
do
    echo "round $i"
    $path/solvers/belos_kosz_mutate.exe $targetpath/testfiles $targetpath/survived $method1 $4 50 --nullfile=ones --tol=1e-12 --max-iters=$iters --belos-xml=$path/solvers/XML/belosdefault.xml --jacobi-xml=$path/solvers/XML/jacobidefault.xml --ilut-xml=$path/solvers/XML/ilutdefault.xml --rhsfile=$path/solvers/RHS/$5.mtx
    rm $targetpath/testfiles/*
    #mkdir $targetpath/survived${i}
    #cp $targetpath/survived/* $targetpath/survived${i}/
    mv $targetpath/survived/*.mtx $targetpath/testfiles/
    
done

