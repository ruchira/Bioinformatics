// clone.h
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
#ifndef CLONE_H
#define CLONE_H
#include "friends_or_foes.h"
#include <map>

using namespace std;

class Clone;

typedef map<Clone *, float> CoefficientEffectOnFitnessMap;

class Clone {
  public:
    Clone() {
      default_fitness = 1.0;
      survival_coefficient = 1.0;
      reproduction_coefficient = 1.0;
    };
    virtual ~Clone() {
    };
    float get_constant_effect_of_clone(Clone *neighbor) const {
      CoefficientEffectOnFitnessMap::const_iterator iter;
      iter = constant_coefficients.find(neighbor);
      if (iter == constant_coefficients.end()) {
        return 0.0;
      } else {
        return iter->second;
      }
    };
    float get_linear_effect_of_clone(Clone *neighbor) const {
      CoefficientEffectOnFitnessMap::const_iterator iter;
      iter = linear_coefficients.find(neighbor);
      if (iter == linear_coefficients.end()) {
        return 0.0;
      } else {
        return iter->second;
      }
    };
    void set_constant_effect_of_clone(Clone *neighbor, 
                                      float constant_coefficient) {
      constant_coefficients[neighbor] = constant_coefficient;
    };
    void set_linear_effect_of_clone(Clone *neighbor, float linear_coefficient) {
      linear_coefficients[neighbor] = linear_coefficient;
    };
    float get_default_fitness(void) const { return default_fitness; };
    void set_default_fitness(float fitness) { default_fitness = fitness; };
    float get_survival_coefficient(void) const { return survival_coefficient; };
    void set_survival_coefficient(float coeff) { 
      survival_coefficient = coeff; 
    };
    float get_reproduction_coefficient(void) const { 
      return reproduction_coefficient;
    };
    void set_reproduction_coefficient(float coeff) {
      reproduction_coefficient = coeff;
    };

  private:
    float default_fitness;
    float survival_coefficient;
    float reproduction_coefficient;
    CoefficientEffectOnFitnessMap constant_coefficients;
    CoefficientEffectOnFitnessMap linear_coefficients;
};

#endif
