import csv
import sys

from os import listdir
from os.path import isfile, join

#read the files located at path (argument of this script) and check if the error in the file (only with 000000 or 111111 BCs) is well equal to the machine epsilon 
#return the number of mistakes

#path 2 your results
path = sys.argv[len(sys.argv)-1]

#absolute tolerance to machine epsilon
tol = 1e-12

#Checking for exactness of results
n_mistake = 0

#list all files in path
files = [f for f in listdir(path) if (isfile(join(path, f)) and ('000000' in f) or ('111111' in f))]


for file in files:
    fcurr = open(path+'/'+file,'r')

    print(path+'/'+file)

    #comparing current results with reference
    for line in csv.reader(fcurr,delimiter=' '):
        buff = list(line)  
        
        err2   = abs(float(buff[1]))
        errInf = abs(float(buff[2]))

        if err2<tol and errInf<tol:
            pass
        else:
            print("WRONG values in "+file+"  ("+buff[1]+" , " + buff[2] + " ... expected 0.0)" )
            n_mistake +=1

    fcurr.close()

print("done, with %d error(s)."%n_mistake)
