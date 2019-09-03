/**
 * @file Profiler.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-08-26
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "Profiler.hpp"

void TimerAgent::start() {
    _count += 1;
    _t0 = MPI_Wtime();
}

void TimerAgent::stop() {
    // get the time
    _t1 = MPI_Wtime();
    // store it
    double dt = _t1 - _t0;
    _timeAcc  = _timeAcc + dt;
    _timeMax  = max(_timeMax, dt);
    _timeMin  = min(_timeMax, dt);
}

void TimerAgent::reset() {
    _t1      = 0.0;
    _t0      = 0.0;
    _timeAcc = 0.0;
    _timeMax = 0.0;
    _timeMin = 0.0;
}

Profiler::Profiler(){
    _name = "default";
}
Profiler::Profiler(string myname){
    _name = myname;
}
Profiler::~Profiler() {
    for (map<string, TimerAgent*>::iterator it = _timeMap.begin(); it != _timeMap.end(); it++) {
        delete (it->second);
    }
}

void Profiler::create(string name) {
    map<string, TimerAgent*>::iterator it = _timeMap.find(name);
    // if it does not already exist
    if (it == _timeMap.end()) {
        _timeMap[name] = new TimerAgent();
        _timeMap[name]->reset();
    }
}

void Profiler::start(string name) {
    _timeMap[name]->start();
}

void Profiler::stop(string name) {
    _timeMap[name]->stop();
}

void Profiler::disp() {
    int commSize, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&commSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

     FILE* file;
    if(rank==0){
        string filename = "prof/"+_name+".prof";
        file = fopen(filename.c_str(),"a+");
        UP_CHECK0(file!=NULL,"Profiler file failed to open");
    }
    
    if(rank ==0){
        printf("=============================================================================================================\n");
        printf("        PROFILER %s  \n",_name.c_str());
        printf("\t-NAME-   \t-Total time-\t-time/call-\t-Min tot time-\t-Max tot time-\t-Mean cnt-\n");

        fprintf(file,"=============================================================================================================\n");
        fprintf(file,"        PROFILER %s  \n",_name.c_str());
        fprintf(file,"\t-NAME-   \t-Total time-\t-time/call-\t-Max tot time-\t-Min tot time-\t-Mean cnt-\n");
    }

    // go through each and compute relevant info
    for (map<string, TimerAgent*>::iterator it = _timeMap.begin(); it != _timeMap.end(); it++) {
        // comnpute counter stuffs
        double localCount =(double) it->second->count();
        double meanCount, maxCount, minCount;
        MPI_Allreduce(&localCount, &meanCount, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localCount, &maxCount, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&localCount, &minCount, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        meanCount /= commSize;

        // comnpute counter stuffs
        double localTime = it->second->timeAcc();
        double meanTime, maxTime, minTime;
        MPI_Allreduce(&localTime, &meanTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&localTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        meanTime /= commSize;

        double meanTimePerCount;
        double localTimePerCount = localTime/localCount;
        MPI_Allreduce(&localTimePerCount,&meanTimePerCount,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        meanTimePerCount /= commSize;

        // if(rank ==0) printf("%15.15s:  \t%02.6f\t%02.6f\t%02.6f\t%09.0f\t%09.0f\t%09.0f\n",it->first.c_str(),meanTime,minTime,maxTime,meanCount,minCount,maxCount);
        if(rank ==0){
            printf("%15.15s:  \t%02.6f\t%02.6f\t%02.6f\t%02.6f\t%09.0f\n",it->first.c_str(),meanTime,meanTimePerCount,minTime,maxTime,meanCount);
            fprintf(file,"%15.15s:  \t%02.6f\t%02.6f\t%02.6f\t%02.6f\t%09.0f\n",it->first.c_str(),meanTime,meanTimePerCount,minTime,maxTime,meanCount);
        }
    }
    if(rank ==0){
        printf("=============================================================================================================\n");
        printf("Total time - the mean of the total time spend in that timer among the processors\n");
        printf("max time - the max total time spend in that timer among the processors\n");
        printf("=============================================================================================================\n");

        fprintf(file,"=============================================================================================================\n");
        fprintf(file,"Total time - the mean of the total time spend in that timer among the processors\n");
        fprintf(file,"max time - the max total time spend in that timer among the processors\n");
        fprintf(file,"=============================================================================================================\n");

        fclose(file);
    }
}