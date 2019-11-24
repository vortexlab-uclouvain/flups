import subprocess
from check_res import check_res

#List of combinations of some boundary conditions in 3 direction:

BCs = [ ["4","4","4","4","4","4"],
        ["4","4","0","4","4","4"],
        ["4","1","1","4","4","4"],
        ["4","1","4","4","4","4"],
        ["4","0","1","4","4","4"],
        ["4","0","1","4","4","1"],
        ["3","3","3","3","3","3"],
        ["4","0","1","4","9","9"],
        ["3","3","3","3","9","9"]]

Kernels = ['0','1','2','3','4']

#Running all combinations of bcs and all kernels
n_success = 0
n_failure = 0

print("Starting the tests...")

i = 0
for bcs in BCs :
    for kern in Kernels :
        i+=1

        print("----- %i -----"%i, flush=True)

        code = ''.join(bcs)

        # if kernel = LGF, we only do the unbounded, if not, we do everything
        # if ((kern=='1' and (bcs==["4","4","4","4","4","4"] or bcs==["3","3","3","3","9","9"])) or (kern != '1') ):
            # Launching test
            #+ ["-oversubscribe"]
        if(bcs[4:6] == ["9","9"]):
            # print("kikouuu from "%i + code)
            r = subprocess.run(["mpirun"] + ["-np"] + ["2"] + ["./flups_validation_nb"] + ["-np"] + ["1"] + ["2"] + ["1"] + ["-k"] + [kern] + ["-res"] + ["16"] + ["16"] + ["1"] + ["-nres"] + ["1"] + ["-bc"] + bcs, capture_output=True)
        else:
            r = subprocess.run(["mpirun"] + ["-np"] + ["2"] + ["./flups_validation_nb"] + ["-np"] + ["1"] + ["2"] + ["1"] + ["-k"] + [kern] + ["-res"] + ["16"] + ["16"] + ["16"] + ["-nres"] + ["1"] + ["-bc"] + bcs, capture_output=True)
        
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
        # else:
        #     print("test %i (BCs : "%i + code + " and k="+ kern+ ") does not apply")

print("%i test succeed out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
