import csv
import os

#read and compare the 'file' located in ./data and ./data_ref. 
#return the number of mistakes
#i is just an index for display
def check_res_3d(i,file):

    tol = 1e-10 #relative error does not exceed tol

    #Checking for exactness of results
    n_mistake = 0

    fcurr = open('./data/'+file,'r')
    
    #creating a dictionnary with reference data
    try: 
        fref  = open('./data_ref/'+file,'r')

        dicref = {}
        for line in csv.reader(fref,delimiter=' '):
            buff = list(line)  
            dicref.update({buff[0] : [float(buff[1]),float(buff[2]),float(buff[3]),float(buff[4]),float(buff[5]),float(buff[6])] })

        fref.close()
    except FileNotFoundError:
        dicref = {}

    #comparing current results with reference
    for line in csv.reader(fcurr,delimiter=' '):
        buff = list(line)  

        vals = dicref.get(buff[0])
        if vals is None:
            n_mistake +=1
            print("test %i: skipped res= "%i +buff[0]+", no ref data.\n     curr: "+buff[1]+" , " + buff[2]+" , "+buff[3]+" , " + buff[4]+" , "+buff[5]+" , " + buff[6] )
            continue
        
        err2   = float('Inf')
        errInf = float('Inf')
        if vals[0] < 1e-14: #absolute error because ref is 0
            err2   = abs((vals[0]-float(buff[1])))
        else:
            err2   = abs((vals[0]-float(buff[1]))/vals[0])

        if vals[1] < 1e-14: #absolute error because ref is 0
            errInf = abs((vals[1]-float(buff[2])))
        else:
            errInf = abs((vals[1]-float(buff[2]))/vals[1])

        if vals[2] < 1e-14: #absolute error because ref is 0
            err2   = err2   + abs((vals[2]-float(buff[3])))
        else:
            err2   = err2   + abs((vals[2]-float(buff[3]))/vals[2])

        if vals[3] < 1e-14: #absolute error because ref is 0    
            errInf = errInf + abs((vals[3]-float(buff[4])))
        else:
            errInf = errInf + abs((vals[3]-float(buff[4]))/vals[3])

        if vals[4] < 1e-14: #absolute error because ref is 0
            err2   = err2   + abs((vals[4]-float(buff[5])))
        else:
            err2   = err2   + abs((vals[4]-float(buff[5]))/vals[4])

        if vals[5] < 1e-14: #absolute error because ref is 0
            errInf = errInf + abs((vals[5]-float(buff[6])))
        else:
            errInf = errInf + abs((vals[5]-float(buff[6]))/vals[5])
            
            
        if err2<(3.0*tol) and errInf<(3.0*tol):
            pass
        else:
            print("test %i: WRONG values for res= "%i +buff[0]+":\n     curr: "+buff[1]+buff[3]+buff[5]+" , " + buff[2]+buff[4]+buff[6]+"\n     ref : %10.12e , %10.12e"%(vals[0],vals[1]))
            n_mistake +=1

    fcurr.close()

    return n_mistake