import subprocess
import csv

#List of combinations of boundary conditions in 1 direction:
BC1s = [["0"],["1"],["4"]]

BC1 = []
for bc1 in BC1s:
    for bc2 in BC1s:
        BC1.append(bc1+bc2)
BC1.append(["3","3"])

Kernels = ['0','1','2','3']

#Running all combinations of bcs
n_success = 0
n_failure = 0

tol = 1e-14

tmp = [["4","4"]]

i = 0
for bcx in BC1 :
    for bcy in BC1 :
        for bcz in BC1:
            i+=1
            code = bcx[0] + bcx[1] + bcy[0] + bcy[1] + bcz[0] + bcz[1]

            #Launching test
            r = subprocess.run(["./flups_validation"] + ["-res"] + ["8"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
            
            if r.returncode != 0 :
                print("test %i (BCs : "%i + code + ") failed with error code ",r.returncode)
                print("=================================== STDOUT =============================================" )
                print(r.stdout.decode())
                print("=================================== STDERR =============================================" )
                print(r.stderr.decode())
                n_failure += 1
                print("=================================== ====== =============================================\n" )
                continue

            #Checking for exactness of results
            n_mistake = 0
            
            fcurr = open('./data/validation_3d_'+code+'_typeGreen=0.txt','r')
            fref  = open('./data_ref/validation_3d_'+code+'_typeGreen=0.txt','r')

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
                    print("test %i (BCs : "%i + code + ") skipped res= "+buff[0]+", no ref data.")
                    continue
                elif abs(vals[0]-float(buff[1]))<tol and abs(vals[1]-float(buff[2]))<tol:
                    pass
                else:
                    n_mistake +=1

            if n_mistake==0:
                print("test %i (BCs : "%i + code + ") succeed")
                n_success += 1    
            else:
                print("test %i (BCs : "%i + code + ") failed with wrong values:")
                n_failure += 1
                fcurr.seek(0)
                fref.seek(0)
                print("=================================== CURRENT VALUES =============================================" )
                print(fcurr.read())
                print("=================================== REFERENCE VALUES ===========================================" )
                print(fref.read())      
                print("=================================== ================ ===========================================\n" )
            
            fcurr.close()
            fref.close()

print("%i test succeed out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
