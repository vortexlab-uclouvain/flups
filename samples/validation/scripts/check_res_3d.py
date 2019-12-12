import csv
import os

#read and compare the 'file' located in ./data and ./data_ref. 
#return the number of mistakes
#i is just an index for display
def check_res_3d(i,file,reffile1,reffile2,reffile3):

    tol = 1e-10 #relative error does not exceed tol

    #Checking for exactness of results
    n_mistake = 0

    fcurr = open('./data/'+file,'r')
    
    #creating a dictionnary with reference data
    try: 
        fref1  = open('./data_ref/'+reffile1,'r')
        fref2  = open('./data_ref/'+reffile2,'r')
        fref3  = open('./data_ref/'+reffile3,'r')

        dicref1 = {}
        dicref2 = {}
        dicref3 = {}
        for line in csv.reader(fref1,delimiter=' '):
            buff = list(line)  
            dicref1.update({buff[0] : [float(buff[1]),float(buff[2])] })
        for line in csv.reader(fref2,delimiter=' '):
            buff = list(line)  
            dicref2.update({buff[0] : [float(buff[1]),float(buff[2])] })
        for line in csv.reader(fref3,delimiter=' '):
            buff = list(line)  
            dicref3.update({buff[0] : [float(buff[1]),float(buff[2])] })

        fref1.close()
        fref2.close()
        fref3.close()
    except FileNotFoundError:
        dicref = {}

    #comparing current results with reference
    for line in csv.reader(fcurr,delimiter=' '):
        buff = list(line)  

        vals1 = dicref.get(buff[0])
        vals2 = dicref.get(buff[0])
        vals3 = dicref.get(buff[0])
        if vals1 is None:
            n_mistake +=1
            print("test %i: skipped res= "%i +buff[0]+", no ref data.\n     curr: "+buff[1]+" , " + buff[2] )
            continue
        if vals2 is None:
            n_mistake +=1
            print("test %i: skipped res= "%i +buff[0]+", no ref data.\n     curr: "+buff[1]+" , " + buff[2] )
            continue
        if vals3 is None:
            n_mistake +=1
            print("test %i: skipped res= "%i +buff[0]+", no ref data.\n     curr: "+buff[1]+" , " + buff[2] )
            continue
        
        err2   = float('Inf')
        errInf = float('Inf')
        if vals1[0] < 1e-14: #absolute error because ref is 0
            err2   = abs((vals1[0]-float(buff[1])))
            errInf = abs((vals2[1]-float(buff[2])))
        if vals2[0] < 1e-14: #absolute error because ref is 0
            err2   = err2   + abs((vals2[0]-float(buff[3])))
            errInf = errInf + abs((vals2[1]-float(buff[4])))
        if vals3[0] < 1e-14: #absolute error because ref is 0
            err2   = err2   + abs((vals3[0]-float(buff[5])))
            errInf = errInf + abs((vals3[1]-float(buff[6])))
        else:
            err2   = abs((vals1[0]-float(buff[1]))/vals1[0])
            errInf = abs((vals1[1]-float(buff[2]))/vals1[1])
            err2   = err2   + abs((vals2[0]-float(buff[3]))/vals2[0])
            errInf = errInf + abs((vals2[1]-float(buff[4]))/vals2[1])
            err2   = err2   + abs((vals3[0]-float(buff[5]))/vals3[0])
            errInf = errInf + abs((vals3[1]-float(buff[6]))/vals3[1])

        if err2<(3.0*tol) and errInf<(3.0*tol):
            pass
        else:
            print("test %i: WRONG values for res= "%i +buff[0]+":\n     curr: "+buff[1]+buff[3]+buff[5]+" , " + buff[2]+buff[4]+buff[6]+"\n     ref : %10.12e , %10.12e"%(vals[0],vals[1]))
            n_mistake +=1

    fcurr.close()

    return n_mistake