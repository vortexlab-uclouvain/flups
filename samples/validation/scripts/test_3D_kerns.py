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

Kernels = ['0','1','2','3','4','5','6','7']

#Running all combinations of bcs and all kernels
n_success = 0
n_failure = 0

print("Starting the tests...")


# centerType = [0, 1]
centerType = [1]
centername = ['NodeCenter', 'CellCenter']
i = 0
for ct in centerType :
    for bcs in BCs :
        for kern in Kernels :
            i+=1

            print("----- %i -----"%i, flush=True)

            code = ''.join(bcs)

            if(kern == '7' and bcs == ["4","0","1","4","9","9"] ) : 
                print("skip kernel 7 and 2dirunbounded... unsupported error due to inherent approximation.")
                continue

            # if kernel = LGF, we only do the unbounded, if not, we do everything
            # if ((kern=='1' and (bcs==["4","4","4","4","4","4"] or bcs==["3","3","3","3","9","9"])) or (kern != '1') ):
                # Launching test
                #+ ["-oversubscribe"]
            if(bcs[4:6] == ["9","9"]):
                r = subprocess.run(["mpirun"] + ["-np"] + ["2"] + ["./flups_validation_nb"] + ["-np"] + ["1"] + ["2"] + ["1"] + ["-k"] + [kern] + ["--center"] + [str(ct)] + ["-res"] + ["16"] + ["16"] + ["1"] + ["-nres"] + ["1"] + ["-bc"] + bcs, capture_output=True)
            else:
                r = subprocess.run(["mpirun"] + ["-np"] + ["2"] + ["./flups_validation_nb"] + ["-np"] + ["1"] + ["2"] + ["1"] + ["-k"] + [kern] + ["--center"] + [str(ct)] + ["-res"] + ["16"] + ["16"] + ["16"] + ["-nres"] + ["1"] + ["-bc"] + bcs, capture_output=True)
            
            if r.returncode != 0 :
                print("test %i ( "%i + centername[ct] + " - BCs : " + code + "with kernel "+kern+") failed with error code ",r.returncode)
                print("=================================== STDOUT =============================================" )
                print(r.stdout.decode())
                print("=================================== STDERR =============================================" )
                print(r.stderr.decode())
                n_failure += 1
                print("=================================== ====== =============================================\n" )
                continue

            #Checking for exactness of results
            print('validation_3d_' + centername[ct] + '_' +code+'_typeGreen='+ kern +'.txt')
            n_mistake = check_res(i,'validation_3d_' + centername[ct] + '_' +code+'_typeGreen='+ kern +'.txt')

            if n_mistake==0:
                print("test %i ( "%i + centername[ct] + " - BCs : " + code + " and k="+ kern+ ") succeeded")
                n_success += 1    
            else:
                print("test %i ( "%i + centername[ct] + " - BCs : " + code + " and k="+ kern+ ") failed with wrong values.")
                print("/!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ -- /!\ \n")
                n_failure += 1
            # else:
            #     print("test %i (BCs : "%i + code + " and k="+ kern+ ") does not apply")

print("%i test succeeded out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
