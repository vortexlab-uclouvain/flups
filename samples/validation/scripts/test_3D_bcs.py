import subprocess
import csv
import sys 
from check_res import check_res

## Check which type of center you try to test
try:
    arg = sys.argv[1]
except IndexError:
    print("/!\ /!\ /!\ WARNING /!\ /!\ /!\ ")
    print("You didn't choose any center type. ")
    print("By default, we will test Cell centered data")
    arg = 1


if(int(arg) == 0):
    print("We will test Node centered data")
    centerType = int(arg)
    centername = 'NodeCenter'
elif(int(arg) == 1): 
    print("We will test Cell centered data")
    centerType = int(arg)
    centername = 'CellCenter'
else: 
    print("/!\ /!\ /!\ WARNING /!\ /!\ /!\ ")
    print("You choose a center type which is not supported. ")
    print("By default, we will test Cell centered data")
    centerType = 1
    centername = 'CellCenter'

#List of combinations of boundary conditions in 1 direction:
BC1s = [["0"],["1"],["4"]]

BC1 = []
for bc1 in BC1s:
    for bc2 in BC1s:
        BC1.append(bc1+bc2)
BC1.append(["3","3"])

#Running all combinations of bcs
n_success = 0
n_failure = 0

BC2 = BC1.copy()
BC2.append(["9","9"])

i = 0
for bcx in BC1 :
    for bcy in BC1 :
        for bcz in BC2:
            i+=1
            code = bcx[0] + bcx[1] + bcy[0] + bcy[1] + bcz[0] + bcz[1]

             #Launching test
            if(bcz == ["9","9"]):
                r = subprocess.run(["./flups_validation_nb"] +["--center"] + [str(centerType)] +["-res"] + ["8"] + ["8"] + ["1"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
            else:
                r = subprocess.run(["./flups_validation_nb"] + ["--center"] + [str(centerType)] +["-res"] + ["8"] + ["8"] + ["8"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
            
            if r.returncode != 0 :
                print("test %i  ( "%i + centername + " - BCs : " + code + ") failed with error code ",r.returncode)
                print("=================================== STDOUT =============================================" )
                print(r.stdout.decode())
                print("=================================== STDERR =============================================" )
                print(r.stderr.decode())
                n_failure += 1
                print("=================================== ====== =============================================\n" )
                continue

             #Checking for exactness of results
            n_mistake = check_res(i,'validation_3d_' + centername + '_' + code +'_typeGreen=0.txt')
            if n_mistake==0:
                print("test %i ( "%i + centername + " - BCs : " + code + ") succeeded")
                n_success += 1    
            else:
                print("test %i ( "%i + centername + " - BCs : " + code + ") failed with wrong values.")
                print("/!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ \n" )
                n_failure += 1

print("%i test succeeded out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
