// /**
//  * @file Profiler.cpp
//  * @author Thomas Gillis and Denis-Gabriel Caprace
//  * @copyright Copyright © UCLouvain 2020
//  * 
//  * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
//  * 
//  * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
//  * 
//  * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
//  * 
//  * Licensed under the Apache License, Version 2.0 (the "License");
//  * you may not use this file except in compliance with the License.
//  * You may obtain a copy of the License at
//  * 
//  *  http://www.apache.org/licenses/LICENSE-2.0
//  * 
//  * Unless required by applicable law or agreed to in writing, software
//  * distributed under the License is distributed on an "AS IS" BASIS,
//  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  * See the License for the specific language governing permissions and
//  * limitations under the License.
//  * 
//  */

// #include "Profiler.hpp"


// /**
//  * @brief Construct a new Timer Agent
//  * 
//  * @param name the name
//  */
// TimerAgent::TimerAgent(string name){
//     name_ = name;
// }

// /**
//  * @brief start the timer
//  * 
//  */
// void TimerAgent::start() {
//     count_ += 1;
//     t0_ = MPI_Wtime();
// }

// /**
//  * @brief stop the timer
//  * 
//  */
// void TimerAgent::stop() {
//     // get the time
//     t1_ = MPI_Wtime();
//     // store it
//     double dt = t1_ - t0_;
//     timeAcc_  = timeAcc_ + dt;
//     timeMax_  = max(timeMax_, dt);
//     timeMin_  = min(timeMax_, dt);
// }

// /**
//  * @brief reset the timer
//  * 
//  */
// void TimerAgent::reset() {
//     t1_      = 0.0;
//     t0_      = 0.0;
//     timeAcc_ = 0.0;
//     timeMax_ = 0.0;
//     timeMin_ = 0.0;
// }

// /**
//  * @brief adds memory to the timer to compute bandwith
//  * 
//  */
// void TimerAgent::addMem(size_t mem){
//     memsize_ += mem;
// }

// /**
//  * @brief add a child to the timer 
//  * 
//  * @param child 
//  */
// void TimerAgent::addChild(TimerAgent* child) {
//     string childName = child->name();
//     map<string, TimerAgent*>::iterator it = children_.find(childName);
//     // if it does not already exist
//     if (it == children_.end()) {
//         children_[childName] = child;
//         child->setDaddy(this);
//     }
// }
// /**
//  * @brief store the dady pointer
//  * 
//  * @param daddy 
//  */
// void TimerAgent::setDaddy(TimerAgent* daddy) {
//     daddy_  = daddy;
//     isroot_ = false;
// }

// /**
//  * @brief display the time accumulated. If it's a ghost timer (no calls), we sum the time of the children
//  * 
//  * @return double 
//  */
// double TimerAgent::timeAcc() const {
//     if (count_ > 0) {
//         return timeAcc_;
//     } else {
//         double sum = 0.0;
//         for (map<string,TimerAgent*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
//             const TimerAgent* child = it->second;
//             sum += child->timeAcc();
//         }
//         return sum;
//     }
// }

// /**
//  * @brief display the min time among all calls. If it's a ghost timer (no calls), we sum the time of the children
//  * 
//  * @return double 
//  */
// double TimerAgent::timeMin() const {
//     if (count_ > 0) {
//         return timeMin_;
//     } else {
//         double sum = 0.0;
//         for (map<string,TimerAgent*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
//             const TimerAgent* child = it->second;
//             sum += child->timeMin();
//         }
//         return sum;
//     }
// }

// /**
//  * @brief display the max time among all calls. If it's a ghost timer (no calls), we sum the time of the children
//  * 
//  * @return double 
//  */
// double TimerAgent::timeMax() const {
//     if (count_ > 0) {
//         return timeMax_;
//     } else {
//         double sum = 0.0;
//         for (map<string,TimerAgent*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
//             const TimerAgent* child = it->second;
//             sum += child->timeMax();
//         }
//         return sum;
//     }
// }

// void TimerAgent::writeParentality(FILE* file, const int level){
//     fprintf(file,"%d;%s",level,name_.c_str());
//     for (map<string, TimerAgent*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
//         const TimerAgent* child = it->second;
//         string childName = child->name();
//         fprintf(file,";%s",childName.c_str());
//     }
//     fprintf(file,"\n");

//     for (map<string, TimerAgent*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
//         it->second->writeParentality(file,level+1);
//     }
// }

// /**
//  * @brief display the time for the TimerAgent
//  * 
//  * @param file 
//  * @param level 
//  * @param totalTime 
//  */
// void TimerAgent::disp(FILE* file,const int level, const double totalTime){
    
//     // check if any proc has called the agent
//     int totalCount;
//     MPI_Allreduce(&count_,&totalCount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
//     // if someone has every call the agent, display it
//     if (totalCount > 0) {
//         // get the size and usefull stuffs
//         int commSize, rank;
//         MPI_Comm_size(MPI_COMM_WORLD, &commSize);
//         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//         double scale = 1.0/commSize;

//         // compute the counters (mean, max, min)
//         double localCount = count_;
//         double meanCount, maxCount, minCount;
//         MPI_Allreduce(&localCount, &meanCount, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         MPI_Allreduce(&localCount, &maxCount, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//         MPI_Allreduce(&localCount, &minCount, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//         meanCount *= scale;

//         // compute times passed inside + children
//         double localTime = timeAcc_;
//         double meanTime, maxTime, minTime;
//         MPI_Allreduce(&localTime, &meanTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         MPI_Allreduce(&localTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//         MPI_Allreduce(&localTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//         meanTime *= scale;
        
//         double meanTimePerCount = meanTime/meanCount;
//         double minTimePerCount = minTime/meanCount;
//         double maxTimePerCount = maxTime/meanCount;
//         double glob_percent = meanTime/totalTime*100.0;

//         // compute the self time  = time passed inside - children
//         double sumChild = 0.0;
//         for (map<string,TimerAgent*>::iterator it = children_.begin(); it != children_.end(); it++) {
//             TimerAgent* child = it->second;
//             sumChild += child->timeAcc();
//         }
//         double locSelfTime = (this->timeAcc()-sumChild);
//         double selfTime;
//         double self_percent;
//         FLUPS_CHECK(locSelfTime >= 0.0,"The timer %s does not include his children",name_.c_str());
//         MPI_Allreduce(&locSelfTime, &selfTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         selfTime *= scale;
//         self_percent = selfTime / totalTime * 100.0;

//         // comnpute the time passed inside the daddy
//         double loc_percent;
//         if (daddy_ != NULL) {
//             double dadLocalTime = daddy_->timeAcc();
//             double dadMeanTime;
//             MPI_Allreduce(&dadLocalTime, &dadMeanTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//             dadMeanTime *= scale;
//             loc_percent = meanTime / dadMeanTime * 100.0;
            
//         } else {
//             loc_percent = 100.0;
//         }

//         // compute the bandwith
//         double localBandTime = timeAcc_;
//         double localBandMemsize = (double) memsize_;
//         double bandMemsize, bandTime, meanBandwidth;
//         MPI_Allreduce(&localBandTime, &bandTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         MPI_Allreduce(&localBandMemsize, &bandMemsize, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//         meanBandwidth =(bandMemsize/bandTime)/std::pow(10.0,6.0);

//         // setup the displayed name
//         string myname = name_;
//         if (level > 1) {
//             myname = " " + myname;
//         }
//         for (int l = 1; l < level; l++) {
//             myname = "--" + myname;
//         }

//         // printf the important information
//         if (rank == 0) {
//             printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", myname.c_str(), glob_percent, loc_percent, meanTime, selfTime, meanTimePerCount, minTimePerCount, maxTimePerCount, meanCount,meanBandwidth);
//             if (file != NULL) {
//                 fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f;%09.2f\n", name_.c_str(), glob_percent, loc_percent, meanTime, selfTime, meanTimePerCount, minTimePerCount, maxTimePerCount, meanCount,meanBandwidth);
//             }
//         }
//     }
//     // recursive call to the childrens
//     for (map<string,TimerAgent*>::iterator it = children_.begin(); it != children_.end(); it++) {
//         TimerAgent* child = it->second;
//         child->disp(file,level+1,totalTime);
//     }
// }

// //===============================================================================================================================
// //===============================================================================================================================
// //===============================================================================================================================


// Profiler::Profiler(): name_("default")
// {
//     createSingle_("root");
// }
// Profiler::Profiler(const string myname): name_(myname)
// {
//     createSingle_("root");
// }
// Profiler::~Profiler() {
//     for (map<string, TimerAgent*>::iterator it = timeMap_.begin(); it != timeMap_.end(); it++) {
//         TimerAgent* current = it->second;
//         delete(current);
//     }
// }

// /**
//  * @brief create a TimerAgent inside the #timeMap_
//  * 
//  * @param name the key of the entry
//  */
// void Profiler::createSingle_(string name) {
//     map<string, TimerAgent*>::iterator it = timeMap_.find(name);
//     // if it does not already exist
//     if (it != timeMap_.end()){
//         timeMap_[name]->reset();
//     }
//     else if (it == timeMap_.end()) {
//         timeMap_[name] = new TimerAgent(name);
//         timeMap_[name]->reset();
//     }
// }

// /**
//  * @brief create a new TimerAgent with "root" as parent
//  * 
//  * @param name the TimerAgent name
//  */
// void Profiler::create(string name) {
//     create(name,"root");
// }

// /**
//  * @brief create a new TimerAgent 
//  * 
//  * @param child the new TimerAgent
//  * @param daddy the dad of the new TimerAgent if it does not exists, it is created
//  */
// void Profiler::create(string child, string daddy) {
//     // create a new guy
//     createSingle_(child);
//     // find the daddy agent in the root
//     map<string, TimerAgent*>::iterator it = timeMap_.find(daddy);
//     if(it == timeMap_.end()){
//         create(daddy);
//     }
//     timeMap_[daddy]->addChild(timeMap_[child]);
// }

// /**
//  * @brief start the timer of the TimerAgent
//  * 
//  * @param name the TimerAgent name
//  */
// void Profiler::start(string name) {
// #ifdef NDEBUG
//     timeMap_[name]->start();
// #else
//     map<string, TimerAgent*>::iterator it = timeMap_.find(name);
//     if (it != timeMap_.end()) {
//         timeMap_[name]->start();
//     }
//     else{
//         // string msg = "timer "+name+ " not found";
//         FLUPS_CHECK(false, "timer %s not found", name.c_str());
//     }
// #endif
// }

// /**
//  * @brief stop the timer of the TimerAgent
//  * 
//  * @param name the TimerAgent name
//  */
// void Profiler::stop(string name) {
// #ifdef NDEBUG
//     timeMap_[name]->stop();
// #else
//     map<string, TimerAgent*>::iterator it = timeMap_.find(name);
//     if (it != timeMap_.end()) {
//         timeMap_[name]->stop();
//     }
//     else{
//         // string msg = "timer "+name+ " not found";
//         FLUPS_CHECK(false, "timer %s not found", name.c_str());
//     }
// #endif
// }

// void Profiler::addMem(string name,size_t mem) {
// #ifdef NDEBUG
//     timeMap_[name]->addMem(mem);
// #else
//     map<string, TimerAgent*>::iterator it = timeMap_.find(name);
//     if (it != timeMap_.end()) {
//         timeMap_[name]->addMem(mem);
//     }
//     else{
//         FLUPS_CHECK(false, "timer %s not found", name.c_str());
//     }
// #endif
// }

// /**
//  * @brief get the accumulated time
//  * 
//  * @param name 
//  * @return double 
//  */
// double Profiler::get_timeAcc(const std::string ref){

//     int commSize;
//     double localTotalTime = timeMap_[ref]->timeAcc();
//     double totalTime;
//     MPI_Comm_size(MPI_COMM_WORLD, &commSize);

//     MPI_Allreduce(&localTotalTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     totalTime /= commSize;

//     return totalTime;
// }


// /**
//  * @brief display the whole profiler using 
//  * 
//  */
// void Profiler::disp() {
//    this->disp("root");
// }
// /**
//  * @brief display the profiler using the timer spent in ref as a reference for the global percentage computation
//  * 
//  * @param ref 
//  */
// void Profiler::disp(const std::string ref) {    
//     int commSize, rank;
//     MPI_Comm_size(MPI_COMM_WORLD, &commSize);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//     //-------------------------------------------------------------------------
//     /** - I/O of the parentality */
//     //-------------------------------------------------------------------------
//     FILE* file;
//     string folder = "./prof";

//     if (rank == 0) {
//         struct stat st = {0};
//         if (stat(folder.c_str(), &st) == -1) {
//                 mkdir(folder.c_str(), 0770);
//         }

//         string filename = folder + "/" + name_ + "_parent.csv";
//         file            = fopen(filename.c_str(), "w+");

//         if (file != NULL) {
//             timeMap_["root"]->writeParentality(file,0);
//             fclose(file);
//         } else {
//             printf("unable to open file %s !", filename.c_str());
//         }
//     }
    

//     //-------------------------------------------------------------------------
//     /** - do the IO of the timing */
//     //-------------------------------------------------------------------------
    
//     if (rank == 0) {
//         string filename = "./prof/" + name_ + "_time.csv";
//         file            = fopen(filename.c_str(), "w+");
//     }
//     // display the header
//     if (rank == 0) {
//         printf("===================================================================================================================================================\n");
//         printf("        PROFILER %s \n", name_.c_str());
//         // printf("\t-NAME-   \t\t\t-%% global-\t-%% local-\t-Total time-\t-Self time-\t-time/call-\t-Min tot time-\t-Max tot time-\t-Mean cnt-\n");
//         printf("%25s|  %-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\n","-NAME-    ", "-% global-", "-% local-", "-Total time-", "-Self time-", "-time/call-", "-Min time-", "-Max time-","-Mean cnt-","-(MB/s)-");
//     }
//     // get the global timing
//     double totalTime = this->get_timeAcc(ref);

//     // display root with the total time
//     timeMap_["root"]->disp(file,0,totalTime);
//     // display footer
//     if (rank == 0) {
//         printf("===================================================================================================================================================\n");
//         printf("%% global - %% of the total time passed inside or in its children (based on the mean time among processors\n");
//         // printf("%% self glob - %% of the total time passed inside (children not included, based on the mean time among processors\n");
//         printf("%% local - %% of the dad's time passed inside or in its children (from the mean time among processors\n");
//         printf("Total time - the total time spend in that timer (averaged among the processors)\n");
//         printf("Self time - the self time spend in that timer = children not included (averaged among the processors)\n");
//         printf("Time/call - the total time spend in that timer per call of the timer (averaged among the processors)\n");
//         printf("Min time - the min time / call spend in that timer among the processors\n");
//         printf("Max time - the max time / call spend in that timer among the processors\n");
//         printf("Mean cnt - the total number of time the timer has been called (averaged among the processors)\n");
//         printf("===================================================================================================================================================\n");
    
//         if (file != NULL) {
//             fclose(file);
//         } else {
//             printf("unable to open file for profiling !");
//         }
//     }
// }
