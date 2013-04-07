// probability.cpp
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
#include "probability.h"
#include "RngStream.h"
#ifdef USING_MPI
#include <mpi.h>
#endif
#include <limits.h>

double clamp_probability(double score) {
  float result = score;
  if (result < 0.0) {
    result = 0.0;
  }
  if (result > 1.0) {
    result = 1.0;
  }
  return result;
}

int Random::num_streams = 1;
int Random::stream_num = 0;
gsl_rng * Random::gsl_rng_ptr = NULL;
RngStream * Random::rngstream_ptr = NULL;
Random * Random::the_random_instance_ptr = NULL;

void setRngStream(void *state, unsigned long int seed) {
  // This ignores the state and just seeds the random number generator from the
  // seeds.
  unsigned long seeds[6];
  for (int i = 0; i < 6; ++i) {
    seeds[i] = seed;
  }
  RngStream::SetPackageSeed(seeds);
}

unsigned long int getRngStream(void *state) {
  if (Random::get_rngstream_ptr() == NULL) {
    throw UninitRngException();
  }
  return Random::get_rngstream_ptr()->RandInt(INT_MIN, INT_MAX);
}

double getDoubleRngStream(void *state) {
  if (Random::get_rngstream_ptr() == NULL) {
    throw UninitRngException();
  }
  return Random::get_rngstream_ptr()->RandU01();
}

static const gsl_rng_type rngstream_type = {
  "RngStream",
  INT_MAX,
  INT_MIN,
  0,  // The RngStream does have an internal state, the Cg array of 6 
      // doubles.  However, there is no need for gsl to mess with this 
      // directly.  So here we say that the size of the state is 0.
  &setRngStream,
  &getRngStream,
  &getDoubleRngStream
};


Random::Random(unsigned long seed) {
  setRngStream(NULL, seed);
  num_streams = 1;
#ifdef USING_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &num_streams);
#endif
  rngstreams = new RngStream[num_streams];
  stream_num = 0;
#ifdef USING_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &stream_num);
#endif
  rngstream_ptr = &(rngstreams[stream_num]);
}

Random::~Random(void) {
  delete[] rngstreams;
  rngstream_ptr = NULL;
  gsl_rng_ptr = NULL;
  num_streams = 0;
  stream_num = -1;
}

void Random::initialize(unsigned long seed) {
  if (Random::the_random_instance_ptr == NULL) {
    the_random_instance_ptr = new Random(seed);
  }
}

void Random::finalize(void) {
  if (Random::the_random_instance_ptr != NULL) {
    delete the_random_instance_ptr;
    the_random_instance_ptr = NULL;
  }
}

bool Random::bernoulli_check(double probability) {
  if (Random::rngstream_ptr == NULL) {
    throw UninitRngException();
  }
  if (rngstream_ptr->RandU01() <= probability) {
    return true;
  } else {
    return false;
  }
}

int Random::rand_int_between_inclusive(int min, int max) {
  if (Random::rngstream_ptr == NULL) {
    throw UninitRngException();
  }
  return rngstream_ptr->RandInt(min, max);
}

double Random::rand_uniform(void) {
  if (Random::rngstream_ptr == NULL) {
    throw UninitRngException();
  }
  return rngstream_ptr->RandU01();
}
