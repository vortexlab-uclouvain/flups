import subprocess
import csv
import sys 
from check_res_3d import check_res_3d

## Check which communication scheme you try to test
try:
    arg = sys.argv[1]
except IndexError:
    print("/!\ /!\ /!\ WARNING /!\ /!\ /!\ ")
    print("You didn't choose any version of the code. ")
    print("By default, we will test the a2a version")
    arg = 'nb'

if(arg != 'isr' and arg!= 'a2a' and arg!='nb'): 
    print("/!\ /!\ /!\ WARNING /!\ /!\ /!\ ")
    print("You choose a version which is not supported. ")
    print("By default, we will test the non blocking version of the code")
    version = 'nb'
else :
    version = arg
    print(f"We will test the {version} version of the code")


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
            print(["./flups_tube_"+version] + ["--np=1,1,1"] + ["--res=16,16,16"] + ["--kernel=0"] + ["--dir="+str(mid[0])] + ["--sym_x="+str(smx[0])] + ["--sym_y="+str(smy[0])])
            r = subprocess.run(["./flups_tube_"+version] + ["--np=1,1,1"] + ["--res=16,16,16"] + ["--kernel=0"] + ["--dir="+str(mid[0])] + ["--sym_x="+str(smx[0])] + ["--sym_y="+str(smy[0])], capture_output=True)

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
            print("checking: vtube_CellCenter_"+code+"_typeGreen=0.txt")
            n_mistake = check_res_3d(i,'vtube_CellCenter_'+code+'_typeGreen=0.txt')

            if n_mistake==0:
                print("test (" + code + ") succeed")
                n_success += 1    
            else:
                print("test (" + code + ") failed with wrong values.")
                print("/!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ \n" )
                n_failure += 1

print("%i test succeed out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
