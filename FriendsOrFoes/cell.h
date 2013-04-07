// cell.h
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
#ifndef CELL_H
#define CELL_H
#include <stdlib.h>
#include "clone.h"

class Cell {
  friend class Population;
  public:
    Cell() : clone_ptr(NULL), alive(false) {};
    virtual ~Cell() {};
    // Accessors
    Clone *get_clone_ptr(void) const { return clone_ptr; };
    virtual float get_volume(void) const = 0;
    bool is_alive(void) const { return alive; };
    float get_previous_fitness(void) const {
      return previous_fitness;
    };
    float get_fitness(void) const { 
      return fitness; 
    };
    void reset_fitness(void) {
      previous_fitness = fitness;
      fitness = clone_ptr->get_default_fitness();
    };
    void increment_fitness_by(float fitness_increment) {
      if (alive) {
        fitness += fitness_increment;
      }
    }
  protected:
    void set_clone_ptr(Clone *a_clone_ptr) { 
      clone_ptr = a_clone_ptr; 
      reset_fitness();
    };
    virtual void envivify(void) { 
      if (clone_ptr != NULL) {
        alive = true; 
        reset_fitness();
      }
    };
    virtual void kill(void) { 
      alive = false; 
      previous_fitness = fitness;
      fitness = 0.0;
    };
    float previous_fitness, fitness;
  private:
    Clone *clone_ptr;
    bool alive;
};

#endif
