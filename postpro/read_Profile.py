import numpy as np
import matplotlib.pyplot as plt

def readProfiler(profName):
    # opent the profiler data
    filename = "../prof/"+ profName + "_time.csv"
    #initialize the time dictionnary
    time = {}
    # Open and read the profiler file
    with open(filename,"r") as file:
        filecont = file.read()
        data = filecont.split("\n",100)
        # read every line
        for s in data:
            entry = s.split(";",100)
            # if the entry has some interesting numbers
            if(len(entry) > 1):
                myname = entry[0]
                # glob_percent, loc_percent,mean_time,self_time,percall_time,min_time,max_time,count
                time[myname] = ([float(entry[1]),float(entry[2]),float(entry[3]),float(entry[4]),float(entry[5]),int(entry[8])])
        # close the file
        file.close()
    return time

def plotProfiler_category(time,profName,catname):
    # get the parentality file
    filename = "../prof/"+ profName + "_parent.csv"
    # create the dictionnary to host the results
    children = list()

    # Open and read the profiler file
    with open(filename,"r") as file:
        filecont = file.read()
        data = filecont.split("\n",100)
        # read every line
        for s in data:
            entry = s.split(";",100)
            # if the entry has some interesting numbers
            length = len(entry)
            if((length> 2) and (entry[1]==catname)):
                for e in entry[2:]:
                    children.append(e)
                #stop reading the file                    
                break
        # close the file
        file.close()
    
    # get the info for each child
    nchild = len(children)
    ind = np.arange(nchild)
    perc = []
    for i in ind:
        child_time = time[children[i]]
        perc.append(child_time[1])

    # do the plot
    plt.barh(ind,perc)
    plt.title(catname)
    plt.yticks(ind,children)
    plt.xlim([0,100])
    plt.grid()
    plt.show()
        

timer = readProfiler("validation")
plotProfiler_category(timer,"validation","solve")