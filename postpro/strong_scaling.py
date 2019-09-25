import numpy as np
import matplotlib.pyplot as plt
import os

#from read_Profile import readProfiler as rp
from read_Profile import *

nx_list = [64]
cpus_list = [1,2,4,8]
thread_list = [1,2,4,8]


waiting = {t: np.empty(0,float) for t in thread_list}
fftw = {t: np.empty(0,float) for t in thread_list}
solve = {t: np.empty(0,float) for t in thread_list}
reorder = {t: np.empty(0,float) for t in thread_list}
cpu = {t: np.empty(0,int) for t in thread_list}
for nx in nx_list:
    for proc in cpus_list:
        for thread in thread_list:
            folder = "../prof"
            name = "validation_res"+str(nx)+"_nrank"+str(proc)+"_nthread"+str(thread)
            # try to read the stuff
            try:
                timer = readProfiler(folder,name)
                print(timer)
                nsolve = timer["solve"][7]

                solve_time = timer["solve"][2]/nsolve
                wait_time = timer["waiting"][2]/nsolve
                fftw_time = timer["fftw"][2]/nsolve
                reorder_time = timer["reorder"][2]/nsolve
                ncore = proc*thread
                solve[thread] = np.append(solve[thread],solve_time)
                waiting[thread] = np.append(waiting[thread],wait_time)
                fftw[thread] = np.append(fftw[thread],fftw_time)
                reorder[thread] = np.append(reorder[thread],reorder_time)
                cpu[thread] = np.append(cpu[thread],ncore)
            except FileNotFoundError:
                continue

# transform the data to np array
# cpu = np.array(cpu)
# cpusec = np.array(cpusec)
# print(time)
# print(cpu)

# the reference is thread=1, first # of cpu
icpuref = 0
ithreadref = 1

for thread in thread_list:
    plt.subplot(2, 2, 1)
    eta = solve[ithreadref][icpuref] / solve[thread]
    plt.plot(cpu[thread], eta, 'o-',label=str(thread)+"shared mem")
    plt.xticks(cpu[thread])

    plt.subplot(2, 2, 2)
    # eta = fftw[thread]
    # eta = (fftw[ithreadref][icpuref]*cpu[ithreadref][icpuref]) / np.multiply(fftw[thread],cpu[thread])
    eta = (solve[ithreadref][icpuref]*cpu[ithreadref][icpuref]) / np.multiply(solve[thread],cpu[thread])
    plt.plot(cpu[thread], eta, 'o-',label=str(thread)+"shared mem")
    plt.xticks(cpu[thread])

    plt.subplot(2, 2, 3)
    # eta = fftw[thread]
    eta = (reorder[ithreadref][icpuref]*cpu[ithreadref][icpuref]) / np.multiply(reorder[thread],cpu[thread])
    plt.plot(cpu[thread], eta, 'o-',label=str(thread)+"shared mem")
    plt.xticks(cpu[thread])

    plt.subplot(2, 2, 4)
    # eta = waiting[thread]
    eta = (waiting[ithreadref][icpuref]*cpu[ithreadref][icpuref]) / np.multiply(waiting[thread],cpu[thread])
    plt.plot(cpu[thread],eta, 'o-',label=str(thread)+"shared mem")
    plt.xticks(cpu[thread])

plt.subplot(2, 2, 1)
plt.plot(cpu[thread], cpu[thread]/cpu[ithreadref][icpuref], 'k--',label="linear")
plt.legend(loc='upper left')
plt.grid()
plt.title("speedup solve")
plt.xlabel("#core")
plt.ylim(0,6)
# plt.xscale('log')
# plt.ylabel("eta strong")
# plt.grid()

plt.subplot(2, 2, 2)
plt.legend()
plt.title("strong scal solve")
plt.xlabel("#core")
plt.xscale('log')
# plt.ylabel("eta strong")
plt.grid()

plt.subplot(2, 2, 3)
plt.legend()
plt.title("strong scal reorder")
plt.xlabel("#core")
plt.xscale('log')
# plt.ylabel("eta strong")
plt.grid()

plt.subplot(2, 2, 4)
plt.legend()
plt.title("strong scal waiting")
plt.xlabel("#core")
plt.xscale('log')
# plt.ylabel("eta strong")
plt.grid()


plt.show()