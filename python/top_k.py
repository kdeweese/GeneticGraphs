#!/usr/bin/python

import numpy
import scipy
import os
import subprocess
import sys

targetpath=sys.argv[1]
k=int(sys.argv[2])
data=numpy.loadtxt(targetpath,usecols=(1,3))
iters=data.T[0]
relres=data.T[1]
data2=numpy.genfromtxt(targetpath,dtype='str',usecols=(0))
filename=data2.T

#inds=numpy.argwhere(iters==maxiters)
#inds=inds.flatten()
#filename=filename[inds]
#relres=(relres[inds])
#iters=iters[inds]

#relres=relres[1:10]
inds=numpy.argsort(iters)
inds=inds.flatten()
relres=relres[inds]
iters=iters[inds]
filename=filename[inds]
for i in xrange(len(iters)-1,len(iters)-k,-1):
    print iters[i]
    print filename[inds]

#print(filename[len(filename)-1])
count=0
for i in xrange(len(filename)-1,len(filename)-k,-1):
#    if(relres[i] < .01):
#        count=count+1

    #subprocess.call(["echo",str(relres[i])])
    basename=os.path.basename(filename[i])
    newfile=sys.argv[3]+basename
    subprocess.call(["cp",filename[i],newfile])

#print(count)
