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

#include <list>
#include "defines.hpp"
#include "mpi.h"

using namespace std;

class TimerAgent {
   protected:
    bool   _isroot  = true;
    int    _count   = 0;
    double _timeAcc = 0.0;
    double _t0      = 0.0;
    double _t1      = 0.0;

    double _timeMax = 0.0;
    double _timeMin = 0.0;

    string _name = "noname";

    TimerAgent*       _daddy = NULL;
    map<string,TimerAgent*> _children;

   public:
    TimerAgent(string name);

    void start();
    void stop();
    void reset();
    void disp(FILE* file,const int level, const double totalTime);

    int    count() const { return _count; };
    bool   isroot() const { return _isroot; };
    string name() const { return _name; };
    double timeAcc() const;
    double timeMin() const;
    double timeMax() const;
    

    void addChild(TimerAgent* child);
    void setDaddy(TimerAgent* daddy);
};

class Profiler {
   protected:
    map<string, TimerAgent*> _timeMap;

    string _name;
    void _createSingle(string name);

   public:
    Profiler();
    Profiler(string myname);
    ~Profiler();

    void create(string name);
    void create(string child, string daddy);

    void start(string name);
    void stop(string name);

    void disp();
};

#endif
