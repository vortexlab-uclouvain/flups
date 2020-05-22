import numpy as np
import matplotlib.pyplot as plt

def readProfiler(folder,profName):
    # opent the profiler data
    filename = folder +"/"+ profName + "_time.csv"
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
                time[myname] = ([float(entry[1]),float(entry[2]),float(entry[3]),float(entry[4]),float(entry[5]),float(entry[6]),float(entry[7]),int(entry[8]),float(entry[9])])
        # close the file
        file.close()
    return time