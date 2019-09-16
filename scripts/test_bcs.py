import subprocess
import filecmp

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

tmp = [["0","0"]]

i = 0
for bcx in BC1 :
    for bcy in tmp :
        for bcz in tmp:
            i+=1
            code = bcx[0] + bcx[1] + bcy[0] + bcy[1] + bcz[0] + bcz[1]
            
            r = subprocess.run(["./flups_validation"] + ["-res"] + ["8"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
            # r = subprocess.run(["mpirun"] + ["-np"] + ["2"] + ["./flups_validation"] + ["-np"] + ["1"] + ["1"] + ["2"] + ["-res"] + ["8"] + ["-bc"] + bcx + bcy + bcz, capture_output=True)
            if r.returncode == 0 :
                f1 = './data/validation_3d_'+code+'_typeGreen=0.err'
                f2 = './data_ref/validation_3d_'+code+'_typeGreen=0.err'
                if(filecmp.cmp(f1,f2)):
                    print("test %i (BCs : "%i + code + ") succeed")
                    n_success += 1    
                else:
                    print("test %i (BCs : "%i + code + ") failed with wrong values:")
                    n_failure += 1
                    print("=================================== CURRENT VALUES =============================================" )
                    print(open(f1,"rb").read().decode('UTF-8'))
                    print("=================================== REFERENCE VALUES =============================================" )
                    print(open(f2,"rb").read().decode('UTF-8'))
            else :
                print("test %i (BCs : "%i + code + ") failed with error code ",r.returncode)
                print("=================================== STDOUT =============================================" )
                print(r.stdout.decode())
                print("=================================== STDERR =============================================" )
                print(r.stderr.decode())
                n_failure += 1
                print("=================================== ====== =============================================\n" )

print("%i test succeed out of %i" % (n_success,n_success+n_failure))
exit(n_failure)
