/**
 * @file Profiler.hpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-08-26
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */
#ifndef PROFILER_HPP
#define PROFILER_HPP
#include <iostream>
#include <map>
#include <string>

#include "mpi.h"

using namespace std;

class TimerAgent {
   protected:
    int    _count   = 0;
    double _timeAcc = 0.0;
    double _t0      = 0.0;
    double _t1      = 0.0;

    double _timeMax = 0.0;
    double _timeMin = 0.0;

   public:
    void start();
    void stop();
    void reset();

    int    count() const{ return _count; };
    double timeAcc() const{ return _timeAcc; };
    double timeMin() const{ return _timeMin; };
    double timeMax() const{ return _timeMax; };
};

class Profiler {
   protected:
    map<string, TimerAgent*> _timeMap;

    string _name;

   public:
   Profiler();
   Profiler(string myname);
    ~Profiler();

    void create(string name);
    void start(string name);
    void stop(string name);

    void disp();
};

#endif
