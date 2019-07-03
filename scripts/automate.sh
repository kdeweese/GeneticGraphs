#!/bin/bash
iterstart=$2
path=/home/kdeweese/geneticbase
targetpath=$1
solver=$3
#mkdir $targetpath/lists
#mkdir $targetpath/histograms
mkdir $targetpath/survived
#$path/scripts/allcg.sh $targetpath $iterstart $solver
#$path/scripts/aggregate.sh $targetpath
#python $path/python/list_process.py $targetpath $iterstart
#rsync -a --delete $path/scripts/empty/ $targetpath/testfiles/
#mv $targetpath/survived/* $targetpath/testfiles/
iters=$iterstart
#for ((iters=$iterstart;iters<=100;iters++))
#do
    #echo "starting iters $iters"
for i in {1..200}
    do
	echo "round $i"
	result=$($path/solvers/belosevolve.exe $targetpath/testfiles $targetpath/survived 5 --nullfile=ones --tol=1e-12 --max-iters=$iters --belos-xml=$path/solvers/XML/belosdefault.xml --precond=$3 --$3-xml=$path/solvers/XML/$3default.xml --rhsfile=$path/solvers/RHS/testrhs.mtx)
	rm $targetpath/testfiles/*
	mv $targetpath/survived/* $targetpath/testfiles/
#	if [ $result -eq 1 ] ;
#	then
#	    let "store = $iters % 5"
#	    if [ $store -eq 0 ]
#	    then
		
#		mkdir $targetpath/survived${iters}
#		cp $targetpath/testfiles/* $targetpath/survived${iters}/
#	    fi
#	    echo "break at iteration $iters round $i"
#	    break
#	fi
	#cp $targetpath/agglist.txt $targetpath/lists/agglist_${i}.txt
	#new=$iters
	#new+=" $i"
	#python $path/python/hist.py $targetpath $new
	#result=$(python $path/python/list_process.py $targetpath $iters)
	#echo "$result"
	
	#if [ $result -eq 0 ] ;
	#then
	#    let "store = $iters % 5"
	#    if [ $store -eq 0 ]
	#    then
	#       mkdir $targetpath/histograms${iters}
	#       mv $targetpath/histograms/* $targetpath/histograms${iters}/
	#       mkdir $targetpath/lists${iters}
	#       mv $targetpath/lists/* $targetpath/lists${iters}/
	#       mkdir $targetpath/survived${iters}
	#       cp $targetpath/survived/* $targetpath/survived${iters}/
	#    else
	#	rm $targetpath/lists/*
	#	rm $targetpath/histograms/*
	#    fi
	#    echo "break at iteration $iters round $i"
	#    break
	#fi
#    done
    #let "store = $iters % 5"
    #if [ $store -eq 0 ]
    #then
    #mkdir $targetpath/histograms${iters}
#	mv $targetpath/histograms/* $targetpath/histograms${iters}/
#	mkdir $targetpath/lists${iters}
#	mv $targetpath/lists/* $targetpath/lists${iters}/
#	mkdir $targetpath/survived${iters}
#	cp $targetpath/survived/* $targetpath/survived${iters}/
 #   else
#	rm $targetpath/lists/*
#	rm $targetpath/histograms/*
 #   fi
  #  echo "continued on iteration $iters round $i"
done
