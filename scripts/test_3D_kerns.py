import subprocess
from check_res import check_res

#List of combinations of some boundary conditions in 3 direction:

BCs = [ ["4","4","4","4","4","4"],
        ["4","4","0","4","4","4"],
        ["4","1","4","4","4","4"],
        ["4","0","1","4","4","4"],
        ["4","0","1","4","4","1"]]

Kernels = ['0','2','3','4']

#Running all combinations of bcs and all kernels
n_success = 0
n_failure = 0

i = 0
for bcs in BCs :
    for kern in Kernels :
        i+=1

        code = ''.join(bcs)

        # Launching test
        r = subprocess.run(["mpirun"] + ["-np"] + ["4"] + ["-oversubscribe"] + ["./flups_validation"] + ["-np"] + ["1"] + ["2"] + ["2"] + ["-k"] + [kern] + ["-res"] + ["16"] + ["-nres"] + ["3"] + ["-bc"] + bcs, capture_output=True)
        
        if r.returncode != 0 :
            print("test %i (BCs : "%i + code + "with kernel "+kern+") failed with error code ",r.returncode)
            print("=================================== STDOUT =============================================" )
            print(r.stdout.decode())
            print("=================================== STDERR =============================================" )
            print(r.stderr.decode())
            n_failure += 1
            print("=================================== ====== =============================================\n" )
            continue

        #Checking for exactness of results
        n_mistake = check_res(i,'validation_3d_'+code+'_typeGreen='+ kern +'.txt')

        if n_mistake==0:
            print("test %i (BCs : "%i + code + " and k="+ kern+ ") succeed")
            n_success += 1    
        else:
            print("test %i (BCs : "%i + code + " and k="+ kern+ ") failed with wrong values.")
            print("/!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ \n")
            n_failure += 1


print("%i test succeed out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
