import scipy
import scipy.io
import sys
import subprocess

filename = sys.argv[1]

A=scipy.io.mmread(filename)
A=A.tocsr()
compos=scipy.sparse.csgraph.connected_components(A,directed=False)
countsize=scipy.zeros(compos[0])
print compos[0]
print compos[1]
if(compos[0]>1):
        print "fail"
        #subprocess.call(["rm",filename])
