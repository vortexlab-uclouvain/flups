import csv
import os

#read and compare the 'file' located in ./data and ./data_ref. 
#return the number of mistakes
#i is just an index for display
def check_res(i, file):

    tol = 1e-14

    #Checking for exactness of results
    n_mistake = 0

    fcurr = open('./data/'+file,'r')
    fref  = open('./data_ref/'+file,'r')

    #creating a dictionnary with reference data
    dicref = {}
    for line in csv.reader(fref,delimiter=' '):
        buff = list(line)  
        dicref.update({buff[0] : [float(buff[1]),float(buff[2])] })

    #comparing current results with reference
    for line in csv.reader(fcurr,delimiter=' '):
        buff = list(line)  
        vals = dicref.get(buff[0])
        if vals is None:
            print("test %i: skipped res= "%i +buff[0]+", no ref data.\n     curr: "+buff[1]+" , " + buff[2] )
            continue
        elif abs(vals[0]-float(buff[1]))<tol and abs(vals[1]-float(buff[2]))<tol:
            pass
        else:
            print("test %i: WRONG values for res= "%i +buff[0]+":\n     curr: "+buff[1]+" , " + buff[2] +"\n     ref : %10.12e , %10.12e"%(vals[0],vals[1]))
            n_mistake +=1

    fcurr.close()
    fref.close()

    return n_mistake