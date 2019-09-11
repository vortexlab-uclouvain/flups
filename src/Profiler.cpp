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


/**
 * @brief Construct a new Timer Agent
 * 
 * @param name the name
 */
TimerAgent::TimerAgent(string name){
    _name = name;
}

/**
 * @brief start the timer
 * 
 */
void TimerAgent::start() {
    _count += 1;
    _t0 = MPI_Wtime();
}

/**
 * @brief stop the timer
 * 
 */
void TimerAgent::stop() {
    // get the time
    _t1 = MPI_Wtime();
    // store it
    double dt = _t1 - _t0;
    _timeAcc  = _timeAcc + dt;
    _timeMax  = max(_timeMax, dt);
    _timeMin  = min(_timeMax, dt);
}

/**
 * @brief reset the timer
 * 
 */
void TimerAgent::reset() {
    _t1      = 0.0;
    _t0      = 0.0;
    _timeAcc = 0.0;
    _timeMax = 0.0;
    _timeMin = 0.0;
}

/**
 * @brief add a child to the timer 
 * 
 * @param child 
 */
void TimerAgent::addChild(TimerAgent* child) {
    string childName = child->name();
    map<string, TimerAgent*>::iterator it = _children.find(childName);
    // if it does not already exist
    if (it == _children.end()) {
        _children[childName] = child;
        child->setDaddy(this);
    }
}
/**
 * @brief store the dady pointer
 * 
 * @param daddy 
 */
void TimerAgent::setDaddy(TimerAgent* daddy) {
    _daddy  = daddy;
    _isroot = false;
}

/**
 * @brief display the time accumulated. If it's a ghost timer (no calls), we sum the time of the children
 * 
 * @return double 
 */
double TimerAgent::timeAcc() const {
    if (_count > 0) {
        return _timeAcc;
    } else {
        double sum = 0.0;
        for (map<string,TimerAgent*>::const_iterator it = _children.begin(); it != _children.end(); it++) {
            const TimerAgent* child = it->second;
            sum += child->timeAcc();
        }
        return sum;
    }
}

/**
 * @brief display the min time among all calls. If it's a ghost timer (no calls), we sum the time of the children
 * 
 * @return double 
 */
double TimerAgent::timeMin() const {
    if (_count > 0) {
        return _timeMin;
    } else {
        double sum = 0.0;
        for (map<string,TimerAgent*>::const_iterator it = _children.begin(); it != _children.end(); it++) {
            const TimerAgent* child = it->second;
            sum += child->timeMin();
        }
        return sum;
    }
}

/**
 * @brief display the max time among all calls. If it's a ghost timer (no calls), we sum the time of the children
 * 
 * @return double 
 */
double TimerAgent::timeMax() const {
    if (_count > 0) {
        return _timeMax;
    } else {
        double sum = 0.0;
        for (map<string,TimerAgent*>::const_iterator it = _children.begin(); it != _children.end(); it++) {
            const TimerAgent* child = it->second;
            sum += child->timeMax();
        }
        return sum;
    }
}

/**
 * @brief display the time for the TimerAgent
 * 
 * @param file 
 * @param level 
 * @param totalTime 
 */
void TimerAgent::disp(FILE* file,const int level, const double totalTime){

    if (_count > 0) {
        // get the size and usefull stuffs
        int commSize, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        double scale = 1.0/commSize;

        // compute the counters (mean, max, min)
        double localCount = _count;
        double meanCount, maxCount, minCount;
        MPI_Allreduce(&localCount, &meanCount, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localCount, &maxCount, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&localCount, &minCount, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        meanCount *= scale;

        // compute time passed inside + children
        double localTime = _timeAcc;
        double meanTime, maxTime, minTime, glob_percent;
        MPI_Allreduce(&localTime, &meanTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&localTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        meanTime *= scale;
        glob_percent = meanTime/totalTime*100.0;

        // compute the self time  = time passed inside - children
        double sumChild = 0.0;
        for (map<string,TimerAgent*>::iterator it = _children.begin(); it != _children.end(); it++) {
            TimerAgent* child = it->second;
            sumChild += child->timeAcc();
        }
        double locSelfTime = (this->timeAcc()-sumChild);
        double selfTime;
        double self_percent;
        FLUPS_CHECK(locSelfTime >= 0.0,"The timer %s does not include his children",_name);
        MPI_Allreduce(&locSelfTime, &selfTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        selfTime *= scale;
        self_percent = selfTime / totalTime * 100.0;

        // compute the time per call
        double meanTimePerCount;
        double localTimePerCount = localTime / localCount;
        MPI_Allreduce(&localTimePerCount, &meanTimePerCount, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        meanTimePerCount *= scale;

        // comnpute the time passed inside the daddy
        double loc_percent;
        if (_daddy != NULL) {
            double dadLocalTime = _daddy->timeAcc();
            double dadMeanTime;
            MPI_Allreduce(&dadLocalTime, &dadMeanTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            dadMeanTime *= scale;
            loc_percent = meanTime / dadMeanTime * 100.0;
            
        } else {
            loc_percent = 100.0;
        }

        // setup the displayed name
        string myname = _name;
        if (level > 1) {
            myname = " " + myname;
        }
        for (int l = 1; l < level; l++) {
            myname = "--" + myname;
        }

        // printf the important information
        if (rank == 0) {
            printf("%-25.25s|  \t%07.4f\t\t%07.4f\t\t%07.4f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%09.0f\n", myname.c_str(), glob_percent,self_percent,loc_percent, meanTime, meanTimePerCount, minTime, maxTime, meanCount);
            fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f\n", myname.c_str(), glob_percent,self_percent,loc_percent, meanTime, meanTimePerCount, minTime, maxTime, meanCount);
        }
    }
    // recursive call to the childrens
    for (map<string,TimerAgent*>::iterator it = _children.begin(); it != _children.end(); it++) {
        TimerAgent* child = it->second;
        child->disp(file,level+1,totalTime);
    }
}

//===============================================================================================================================
//===============================================================================================================================
//===============================================================================================================================


Profiler::Profiler(){
    _name = "default";
    _createSingle("root");
}
Profiler::Profiler(string myname){
    _name = myname;
    _createSingle("root");
}
Profiler::~Profiler() {
    for (map<string, TimerAgent*>::iterator it = _timeMap.begin(); it != _timeMap.end(); it++) {
        delete (it->second);
    }
}

/**
 * @brief create a TimerAgent inside the #_timeMap
 * 
 * @param name the key of the entry
 */
void Profiler::_createSingle(string name) {
    map<string, TimerAgent*>::iterator it = _timeMap.find(name);
    // if it does not already exist
    if (it != _timeMap.end()){
        _timeMap[name]->reset();
    }
    else if (it == _timeMap.end()) {
        _timeMap[name] = new TimerAgent(name);
        _timeMap[name]->reset();
    }
}

/**
 * @brief create a new TimerAgent with "root" as parent
 * 
 * @param name the TimerAgent name
 */
void Profiler::create(string name) {
    create(name,"root");
}

/**
 * @brief create a new TimerAgent 
 * 
 * @param child the new TimerAgent
 * @param daddy the dad of the new TimerAgent if it does not exists, it is created
 */
void Profiler::create(string child, string daddy) {
    // create a new guy
    _createSingle(child);
    // find the daddy agent in the root
    map<string, TimerAgent*>::iterator it = _timeMap.find(daddy);
    if(it == _timeMap.end()){
        create(daddy);
    }
    _timeMap[daddy]->addChild(_timeMap[child]);
}

/**
 * @brief start the timer of the TimerAgent
 * 
 * @param name the TimerAgent name
 */
void Profiler::start(string name) {
#ifdef NDEBUG
    _timeMap[name]->start();
#else
    map<string, TimerAgent*>::iterator it = _timeMap.find(name);
    if (it != _timeMap.end()) {
        _timeMap[name]->start();
    }
    else{
        string msg = "timer "+name+ " not found";
        FLUPS_ERROR(msg);
    }
#endif
}

/**
 * @brief stop the timer of the TimerAgent
 * 
 * @param name the TimerAgent name
 */
void Profiler::stop(string name) {
#ifdef NDEBUG
    _timeMap[name]->stop();
#else
    map<string, TimerAgent*>::iterator it = _timeMap.find(name);
    if (it != _timeMap.end()) {
        _timeMap[name]->stop();
    }
    else{
        string msg = "timer "+name+ " not found";
        FLUPS_ERROR(msg);
    }
#endif
}

/**
 * @brief display the whole profiler
 * 
 */
void Profiler::disp() {
    int commSize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE* file;
    if (rank == 0) {
        string filename = "prof/" + _name + ".csv";
        file            = fopen(filename.c_str(), "w+");
    }
    // display the header
    if (rank == 0) {
        printf("===================================================================================================================================================\n");
        printf("        PROFILER %s  \n", _name.c_str());
        printf("\t-NAME-   \t\t\t-%% global-\t-self %% glob-\t-%% local-\t-Total time-\t-time/call-\t-Min tot time-\t-Max tot time-\t-Mean cnt-\n");
    }
    // get the global timing
    double localTotalTime = _timeMap["root"]->timeAcc();
    double totalTime;
    MPI_Allreduce(&localTotalTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    totalTime /= commSize;
    // display root with the total time
    _timeMap["root"]->disp(file,0,totalTime);
    // display footer
    if (rank == 0) {
        printf("===================================================================================================================================================\n");
        printf("%% global - %% of the total time passed inside or in its children (based on the mean time among processors\n");
        printf("%% self glob - %% of the total time passed inside (children not included, based on the mean time among processors\n");
        printf("%% local - %% of the dad's time passed inside or in its children (from the mean time among processors\n");
        printf("Total time - the total time spend in that timer (averaged among the processors)\n");
        printf("Time/call - the total time spend in that timer per call of the timer (averaged among the processors)\n");
        printf("Min time - the min total time spend in that timer among the processors\n");
        printf("Max time - the max total time spend in that timer among the processors\n");
        printf("Mean cnt - the total number of time the timer has been called (averaged among the processors)\n");
        printf("===================================================================================================================================================\n");
        fclose(file);
    }
}