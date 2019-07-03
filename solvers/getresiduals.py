import string
import sys



f1 = open(sys.argv[1],'r')
f2 = open('residuals.txt','w')


count = 0

while True:
    line = f1.readline()
    if line == "":
        break


    index1 = string.find(line, "residual [ 0 ] = ")
    index2 = string.find(line, "<")
    

    if index1 != -1:
        if index2 == -1:
            index2 = string.find(line, ">")

        #substring = line[temp+17]
        f2.write(str(count))
        f2.write(" ")
        f2.write(line[index1+17:index2]+'\n')
    
        count = count + 1
        
        

f1.close()
f2.close()

