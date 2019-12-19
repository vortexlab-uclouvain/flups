import subprocess
import csv
from check_res_3d import check_res_3d

#List of combinations of boundary conditions in 1 direction:
dirv = [["0"],["1"],["2"]]
sym = [["0"],["1"],["-1"]]

#Running all combinations of bcs
n_success = 0
n_failure = 0

print("Starting the tests...")

i = 0
for mid in dirv :
    for smx in sym :
        for smy in sym:
            i+=1
            code=mid[0]+smx[0]+smy[0]

            print("----- %i -----"%i, flush=True)
            #Launching test
            r = subprocess.run(["./flups_tube_a2a"] + ["-np"] + ["1"] + ["1"] + ["1"] + ["-res"] + ["16"] + ["16"] + ["16"] + ["-k"] + ["0"] + ["-d"] + mid + ["-smx"] + smx + ["-smy"] + smy, capture_output=True)
            
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
            print("checking: vtube_"+code+"_typeGreen=0.txt")
            n_mistake = check_res_3d(i,'vtube_'+code+'_typeGreen=0.txt')

            if n_mistake==0:
                print("test (" + code + ") succeed")
                n_success += 1    
            else:
                print("test (" + code + ") failed with wrong values.")
                print("/!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ \n" )
                n_failure += 1

print("%i test succeed out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
