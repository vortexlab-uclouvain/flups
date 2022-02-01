import subprocess
import csv
from check_res import check_res

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

centerType = [0, 1]
centername = ['NodeCenter', 'CellCenter']

i = 0
for ct in centerType :
    for bcx in BC1 :
        for bcy in BC1 :
            for bcz in BC2:
                i+=1
                code = bcx[0] + bcx[1] + bcy[0] + bcy[1] + bcz[0] + bcz[1]

                #Launching test
                if(bcz == ["9","9"]):
                    r = subprocess.run(["./flups_validation_nb"] +["--center"] + [str(ct)] +["-res"] + ["8"] + ["8"] + ["1"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
                else:
                    r = subprocess.run(["./flups_validation_nb"] + ["--center"] + [str(ct)] +["-res"] + ["8"] + ["8"] + ["8"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
                
                if r.returncode != 0 :
                    print("test %i  ( "%i + centername[ct] + " - BCs : " + code + ") failed with error code ",r.returncode)
                    print("=================================== STDOUT =============================================" )
                    print(r.stdout.decode())
                    print("=================================== STDERR =============================================" )
                    print(r.stderr.decode())
                    n_failure += 1
                    print("=================================== ====== =============================================\n" )
                    continue

                #Checking for exactness of results
                n_mistake = check_res(i,'validation_3d_' + centername[ct] + '_' + code +'_typeGreen=0.txt')
                if n_mistake==0:
                    print("test %i ( "%i + centername[ct] + " - BCs : " + code + ") succeeded")
                    n_success += 1    
                else:
                    print("test %i ( "%i + centername[ct] + " - BCs : " + code + ") failed with wrong values.")
                    print("/!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ \n" )
                    n_failure += 1

print("%i test succeeded out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
