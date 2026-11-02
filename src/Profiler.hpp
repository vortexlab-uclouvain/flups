/**
 * @file Profiler.hpp
 * @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

// #ifndef PROFILER_HPP
// #define PROFILER_HPP
// #include <iostream>
// #include <map>
// #include <string>
// #include <sys/types.h>
// #include <sys/stat.h>
// #include <unistd.h>

// #include <list>
// #include "defines.hpp"
// #include "mpi.h"

// using namespace std;

// #if defined(PROF)
// #define PROF_START(name) if (prof_ != NULL) prof_->start(name);
// #define PROF_STARTi(name,ip) if (prof_ != NULL) prof_->start(name+to_string(ip));
// #define PROF_STOP(name) if (prof_ != NULL) prof_->stop(name);
// #define PROF_STOPi(name,ip) if (prof_ != NULL) prof_->stop(name+to_string(ip));
// #else
// #define PROF_START(name) ((void)0);
// #define PROF_STARTi(name,ip) ((void)0);
// #define PROF_STOP(name) ((void)0);
// #define PROF_STOPi(name,ip) ((void)0);
// #endif

// class TimerAgent {
//    protected:
//     bool   isroot_  = true;
//     int    count_   = 0;
//     double timeAcc_ = 0.0;
//     double t0_      = 0.0;
//     double t1_      = 0.0;
//     size_t memsize_ = 0;

//     double timeMax_ = 0.0;
//     double timeMin_ = 0.0;

//     string name_ = "noname";

//     TimerAgent*       daddy_ = NULL;
//     map<string,TimerAgent*> children_;

//    public:
//     TimerAgent(string name);

//     void start();
//     void stop();
//     void reset();
//     void addMem(size_t mem);
//     void disp(FILE* file, const int level, const double totalTime);

//     int    count() const { return count_; };
//     bool   isroot() const { return isroot_; };
//     string name() const { return name_; };
//     double timeAcc() const;
//     double timeMin() const;
//     double timeMax() const;

//     void addChild(TimerAgent* child);
//     void setDaddy(TimerAgent* daddy);
//     void writeParentality(FILE* file, const int level);
// };

// class Profiler {
//    protected:
//     map<string, TimerAgent*> timeMap_;

//     const string name_;
//     void createSingle_(string name);

//    public:
//     Profiler();
//     Profiler(const string myname);
//     ~Profiler();

//     void create(string name);
//     void create(string child, string daddy);

//     void start(string name);
//     void stop(string name);
//     void addMem(string name,size_t mem);

//     double get_timeAcc(const std::string ref);

//     void disp();
//     void disp(const std::string ref);
// };

// #endif
