// probability.h
// Author: Ruchira S. Datta
// Copyright (c) 2013, Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// o Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// o Redistributions in binary form mus reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// o Neither the name of the University of California, San Francisco nor the
// names of its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PUPROSE ARE
// DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMTIED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
#ifndef PROBABILITY_H
#define PROBABILITY_H
#include <gsl/gsl_rng.h>
#include <string>
#include <exception>
#include "RngStream.h"

class UninitRngException : public std::exception {
  public:
    UninitRngException(std::string m
                    ="Tried to use uninitialized random number generator.\n") :
                   msg(m) {};
    virtual ~UninitRngException() throw() {};
    const char *what() const throw() { return msg.c_str(); };
  private:
    std::string msg;
};

class Random {
  public:
    static void initialize(unsigned long seed = 12345);
    static void finalize(void);
    static int rand_int_between_inclusive(int min, int max);
    static double rand_uniform();
    static bool bernoulli_check(double probability);
    static unsigned int binomial_check(double probability, unsigned int n);
    static gsl_rng *get_gsl_rng_ptr(void) { return gsl_rng_ptr; }
    static RngStream *get_rngstream_ptr(void) { return rngstream_ptr; }
  private:
    Random(unsigned long seed=12345);
    virtual ~Random();
    static Random *the_random_instance_ptr;
    static int num_streams;
    static int stream_num;
    static RngStream *rngstream_ptr;
    static gsl_rng *gsl_rng_ptr;
    RngStream *rngstreams;
};

extern double clamp_probability(double score);
#endif
