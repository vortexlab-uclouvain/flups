/**
 * @file Profiler.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2020
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef PROFILER_HPP
#define PROFILER_HPP
#include <iostream>
#include <map>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <list>
#include "defines.hpp"
#include "mpi.h"

using namespace std;

#if defined(PROF)
#define PROF_START(name) if (prof_ != NULL) prof_->start(name);
#define PROF_STARTi(name,ip) if (prof_ != NULL) prof_->start(name+to_string(ip));
#define PROF_STOP(name) if (prof_ != NULL) prof_->stop(name);
#define PROF_STOPi(name,ip) if (prof_ != NULL) prof_->stop(name+to_string(ip));
#else
#define PROF_START(name) ((void)0);
#define PROF_STARTi(name,ip) ((void)0);
#define PROF_STOP(name) ((void)0);
#define PROF_STOPi(name,ip) ((void)0);
#endif

class TimerAgent {
   protected:
    bool   isroot_  = true;
    int    count_   = 0;
    double timeAcc_ = 0.0;
    double t0_      = 0.0;
    double t1_      = 0.0;
    size_t memsize_ = 0;

    double timeMax_ = 0.0;
    double timeMin_ = 0.0;

    string name_ = "noname";

    TimerAgent*       daddy_ = NULL;
    map<string,TimerAgent*> children_;

   public:
    TimerAgent(string name);

    void start();
    void stop();
    void reset();
    void addMem(size_t mem);
    void disp(FILE* file, const int level, const double totalTime);

    int    count() const { return count_; };
    bool   isroot() const { return isroot_; };
    string name() const { return name_; };
    double timeAcc() const;
    double timeMin() const;
    double timeMax() const;

    void addChild(TimerAgent* child);
    void setDaddy(TimerAgent* daddy);
    void writeParentality(FILE* file, const int level);
};

class Profiler {
   protected:
    map<string, TimerAgent*> timeMap_;

    const string name_;
    void createSingle_(string name);

   public:
    Profiler();
    Profiler(const string myname);
    ~Profiler();

    void create(string name);
    void create(string child, string daddy);

    void start(string name);
    void stop(string name);
    void addMem(string name,size_t mem);

    double get_timeAcc(const std::string ref);

    void disp();
    void disp(const std::string ref);
};

#endif
