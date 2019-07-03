#!/usr/bin/python

import numpy
import scipy
import os
import subprocess
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg')
maxiters=int(sys.argv[1])
histname=sys.argv[2]
print type(maxiters)
data=numpy.loadtxt("lists/agglist_"+histname+".txt",usecols=(1,3))
iters=data.T[0]
relres=data.T[1]
data2=numpy.genfromtxt("lists/agglist_"+histname+".txt",dtype='str',usecols=(0))
filename=data2.T
inds=numpy.argwhere(iters==maxiters)
inds=inds.flatten()
filename=filename[inds]
relres=(relres[inds])
iters=iters[inds]

#print relres
print relres.max()
print relres.min()
#numpy.histogram(relres)
plt.figure(figsize=(14,5))
plt.xlabel('Residual Error Bins')
plt.ylabel('Count')
plt.title("Histogram of Generation " + histname + " Fitness (Population " + str(len(iters)) + ")")
list=[]
start=-12
for i in xrange(0,151):
    list.append(start)
    start=start+.1
alist=numpy.array(list)
binlist=10**alist
#binlist=[0,1e-14,.316e-13,1e-13,.316e-12,1e-12,.316e-11,1e-11,.316e-10,1e-10,.316e-9,1e-9,.316e-8,1e-8,.316e-7,1e-7,.316e-6,1e-6,.316e-5,1e-5,.316e-4,1e-4,.316e-3,1e-3,.316e-2,1e-2,.316e-1,1e-1,.316e0,1e0,.316e1,1e1]
plt.hist(relres,bins=binlist,log=True,range=(0,1))
plt.xscale('log')
plt.ylim(0.1,.5e5)
plt.xlim(1e-12,1.1)
#print n
#print bins
#print patches
#plt.show()
plt.savefig("histograms/hist_"+histname+".png")
