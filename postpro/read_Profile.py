import numpy as np
import matplotlib.pyplot as plt

def readProfiler(profName):
    # opent the profiler data
    filename = "../prof/"+ profName + ".csv"
    # create the dictionnary to host the results
    keys = ["name","tot","percall","count"]
    time = {k: [] for k in keys}

    # Open and read the profiler file
    with open(filename,"r") as file:
        filecont = file.read()
        data = filecont.split("\n",100)
        # read every line
        for s in data:
            entry = s.split(";",100)
            # if the entry has some interesting numbers
            if(len(entry) > 1):
                time.get("name").append(entry[0])
                time.get("tot").append(float(entry[1]))
                time.get("percall").append(float(entry[2]))
                time.get("count").append(int(entry[5]))
        # close the file
        file.close()
    return time


def plotProfiler_row(time,name):
    func        = time.get("name")
    meantime    = time.get("tot")
    # timepercall = time.get("percall")

    ind = np.arange(len(func))
    p1  = plt.barh(ind,meantime)

    plt.title(name)
    plt.yticks(ind,func)
    plt.show()

def plotProfiler_percent(time,name):
    func        = time.get("name")
    totTime    = np.array(time.get("tot"))
    #get the different categories
    category = ["solve","reorder"]

    
    cat_name    = {k: [] for k in category}
    cat_totTime = {k: [] for k in category}
    
    # get the colormap
    cmap = plt.get_cmap("tab10")
    ind = range(len(category))
    # for each category
    for ic in ind:
        cat = category[ic]
        start = 0
        # search among the function names
        for i in range(len(func)):
            # if the function name matches, add the time and the name
            if cat+"_" in func[i]:
                cat_name[cat].append(func[i])
                cat_totTime[cat].append(totTime[i])

        # convert to arrays
        cat_totTime[cat] = np.array(cat_totTime[cat])
        # compute the cumulative sum
        mycum = cat_totTime[cat].cumsum()
        mymax = mycum[-1]
        
        # do the percent plot
        for i in range(len(cat_name[cat])):
            width = cat_totTime[cat][i]/mymax
            left = start
            plt.barh(ic,width,left=left,color=cmap(i))
            start = width
    
    plt.yticks(ind,category)
    plt.show()


timer = readProfiler("FFTW_Solver_128_128_128")
# plotProfiler_row(timer,"FFTW_Solver_128_128_128")
plotProfiler_percent(timer,"FFTW_Solver_128_128_128")