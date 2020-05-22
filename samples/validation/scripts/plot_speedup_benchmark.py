import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as tck
import matplotlib.lines as lines
import tikzplotlib
import os
import sys

from read_Profile import *

#########################################################################
# Create weak scalability plots using output of FLUPS profiler
#
#########################################################################

#path 2 your results (*time.csv)
path = sys.argv[len(sys.argv)-1]

nx_list   = [128] #resolution (should be a 1-entry table)
bcs       = ["000000"] #boundary conditions types, among ["000000"]|["111111"]|["144044"]
cpus_list = [256,512,1024] #list of comm sizes

savefig   = 0 #0|1, save figs as PNGs

#########################################################################

# Fixed parameters

comm_list = ["a2a"] #["a2a","nb"]
thread_list = [1]

linestyle_v = {"a2a": "o-", "nb": "o-"}
mycolor = plt.get_cmap('tab20c')
colorstyle_v   = {"a2a":{"000000":mycolor(0),"111111": mycolor(1),"144044": mycolor(2)},"nb":{"000000":mycolor(0),"111111": mycolor(1),"144044": mycolor(2)}}
colorstyle_vth = {"a2a":{1: mycolor(0),4: mycolor(1)},"nb":{1: mycolor(0),4: mycolor(1)}}

min_core = 1
# the reference is thread=1, first # of cpu
ithreadref = 1 # we take one thraead as ref
icpuref = 0
comm_ref = "a2a"
bc_ref = "none"

type1 = "speedup" #"speedup"|"eta"

#init
waiting = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
fftw = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
solve = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
reorder = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
mem2buf = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
buf2mem = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
cpu = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}
bandwidth = {t: {t: {t: {t: np.empty(0,float) for t in thread_list} for t in comm_list} for t in nx_list} for t in bcs}

for bc in bcs:
    for comm in comm_list:
        for nx in nx_list:
            for procx in cpus_list:
                for thread in thread_list:
                    
                    folder = path
                    name = "prof_bc"+ bc +"_res"+str(nx)+"_nrank"+str(int(procx/thread))+"_nthread"+str(thread)

                    # try to read the stuff
                    try:
                        timer = readProfiler(folder,name)
                        print("reading " + folder + name)
                        nsolve = timer["solve"][7]

                        solve_time = timer["solve"][2]/nsolve

                        wait_time = 0.0
                        my_band = 0.0
                        for k in timer.keys():
                            if "all_2_all" in k and not ("all_2_all_v" in k):
                                wait_time = wait_time+timer[k][2]/nsolve
                                my_band = my_band+timer[k][8]
                                print(k + " " + str(wait_time))
                            if "all_2_all_v" in k:
                                wait_time = wait_time+timer[k][2]/nsolve
                                my_band = my_band+ timer[k][8]
                                print(k + " " + str(wait_time))
                            if "waiting" in k:
                                wait_time = wait_time+timer[k][2]/nsolve
                                my_band = my_band+ timer[k][8]
                                print(k + " " + str(wait_time))
                            if "mem2buf" in k:
                                m2b_time = timer[k][2]/nsolve
                            if "buf2mem" in k:
                                b2m_time = m2b_time = timer[k][2]/nsolve  

                        fftw_time = timer["fftw"][2]/nsolve
                        reorder_time = timer["reorder"][2]/nsolve
                        

                        ncore = procx #*thread
                        if(ncore >= min_core):                            
                            solve[bc][nx][comm][thread] = np.append(solve[bc][nx][comm][thread],solve_time)
                            waiting[bc][nx][comm][thread] = np.append(waiting[bc][nx][comm][thread],wait_time)
                            fftw[bc][nx][comm][thread] = np.append(fftw[bc][nx][comm][thread],fftw_time)
                            reorder[bc][nx][comm][thread] = np.append(reorder[bc][nx][comm][thread],reorder_time)
                            mem2buf[bc][nx][comm][thread] = np.append(mem2buf[bc][nx][comm][thread],m2b_time)
                            buf2mem[bc][nx][comm][thread] = np.append(buf2mem[bc][nx][comm][thread],b2m_time)
                            cpu[bc][nx][comm][thread] = np.append(cpu[bc][nx][comm][thread],ncore)
                            bandwidth[bc][nx][comm][thread] = np.append(bandwidth[bc][nx][comm][thread],my_band)
                    except FileNotFoundError:
                        print("could not find profiler file for "+folder+name)
                        continue


myticks= np.empty(0,int)
coreticks= np.empty(0,int)

barcol = plt.get_cmap('tab20b')
barcolc = plt.get_cmap('tab20c')

for nx in nx_list:

    fig0 = plt.figure(constrained_layout=True)
    ax0 = plt.gca()
    fig2 = plt.figure(constrained_layout=True)
    ax2 = plt.gca()
    
    for bc in bcs:

        for comm in comm_list:
            linestyle = linestyle_v[comm]


            #==================================================
            for thread in thread_list:

                if len(thread_list)>1:
                    # legd= "shm"+str(thread)+" - "+comm
                    legd= "shm"+str(thread)
                else:
                    legd=bc+" - "+comm

                #----------------------
                # PLOT 1
                if type1=="eta":
                    if comm_ref=="none":
                        c_ref = comm
                    else:
                        c_ref = comm_ref
                    if ithreadref=="none":
                        t_ref=thread
                    else:
                        t_ref=ithreadref
                    if bc_ref=="none":
                        m_ref=bc
                    else:
                        m_ref=bc_ref

                    eta = solve[m_ref][nx][c_ref][t_ref][icpuref]/solve[bc][nx][comm][thread]
                else:
                    eta = solve[bc][nx][comm][thread]

                if len(thread_list)>1:
                    ax0.plot(cpu[bc][nx][comm][thread], eta, linestyle,label=legd,color=colorstyle_vth[comm][thread])
                else:
                    ax0.plot(cpu[bc][nx][comm][thread], eta, linestyle,label=legd,color=colorstyle_v[comm][bc])
                ax0.set_xscale('log')
                ax0.minorticks_off()

                myticks = np.union1d(myticks,cpu[bc][nx][comm][thread])
                ax0.set_xticks(myticks)
                ax0.set_xticklabels(myticks.astype(int))

                #----------------------
                # PLOT 3
                width = 0.9 # width of the bars
                # plot for each CPU + a shift
                
                cores = cpu[bc][nx][comm][thread]
                x = (np.log(cores)-np.log(min_core))/np.log(2) - 1
                
                tmp_fftw = fftw[bc][nx][comm][thread]
                tmp_waiting = waiting[bc][nx][comm][thread]
                tmp_mem2buf = mem2buf[bc][nx][comm][thread]
                tmp_buf2mem = buf2mem[bc][nx][comm][thread]
                # tmp_reorder = reorder[nx][comm][thread] - tmp_waiting - tmp_mem2buf - tmp_buf2mem
                tmp_solve = solve[bc][nx][comm][thread]
                tmp_misc = tmp_solve - tmp_fftw - tmp_waiting - tmp_mem2buf - tmp_buf2mem

                # plot the shit
                bh6=ax2.bar(x,tmp_misc,width=width,bottom=tmp_fftw+tmp_waiting+tmp_mem2buf+tmp_buf2mem,color=barcolc(17),alpha=1,edgecolor='w',label="misc")
                bh4=ax2.bar(x,tmp_buf2mem,width=width,bottom=tmp_fftw+tmp_waiting+tmp_mem2buf,color=barcol(11),alpha=1,edgecolor='w',label="buf to mem")
                bh3=ax2.bar(x,tmp_mem2buf,width=width,bottom=tmp_fftw+tmp_waiting,color=barcol(10),alpha=1,edgecolor='w',label="mem to buf")
                bh2=ax2.bar(x,tmp_waiting,width=width,bottom=tmp_fftw,color=barcol(9),alpha=1,edgecolor='w',label="MPI comm")
                bh1=ax2.bar(x,tmp_fftw,width=width,color=barcol(4),alpha=1,edgecolor='w',label="fftw")
                

                # set the ticks and the limits
                myticks = np.union1d(myticks,cores)
                coreticks = np.union1d(coreticks,x)
                ax2.set_xticks(coreticks)
                ax2.set_xticklabels(myticks.astype(int))
                # ax2.set_yticks(tmp_solve)
                # ax2.set_yticklabels(tmp_solve.astype(float))
                # ax2.set_ylim(0,3.5)
                ax2.set_xlabel("#core")
                ax2.set_ylabel("[s]")
                # plt.show()

                # do the io of the bar graph
                ax2.yaxis.grid(b=True,which='major')
                ax2.yaxis.grid(b=True,which='minor', linestyle='--')
                # ax2.legend((bh1,bh2,bh3,bh4,bh6),("FFT","COMM-MPI","COMM-MEM2BUF","COMM-BUF2MEM","MISC"),loc='lower right')
                ax2.legend(loc='lower right')
                
                if savefig==1:
                    plt.savefig(path+"/weak_timepersolve_"+ str(nx) +"_"+comm+"_shm"+str(thread)+"_"+bc+".png",format="png")
                
    #==================================================
    #-------- PLOT 0,0

    ax0.legend(loc='lower left')
    # ax0.set_title("[SOLVE] weak scaling")
    ax0.set_xlabel("N_{cpu}")
    if type1=="eta":
        ax0.set_ylabel("$\eta_{weak}$")
        fname = "scaling"
    else:
        ax0.set_ylabel("$t_{weak} [s]$")
        fname = "speedup"
    ax0.set_ylim(0,)
    ax0.yaxis.grid()
    
    if savefig==1:
        plt.figure(1)
        plt.savefig(path+"/weak_"+fname+"_"+ str(nx) +".png",format="png")
    
    plt.show()
