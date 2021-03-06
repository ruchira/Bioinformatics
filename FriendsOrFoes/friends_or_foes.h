// friends_or_foes.h
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
#ifndef FRIENDS_OR_FOES_H
#define FRIENDS_OR_FOES_H

#include "Poco/Util/Application.h"
#include "Poco/Util/Option.h"
#include "Poco/Util/OptionSet.h"
#include "clone.h"
#include "population.h"
#include "hex_population.h"

#include <iostream>
#include <vector>
using namespace Poco::Util;

class FriendsOrFoesApp: public Application {
public:
  FriendsOrFoesApp(): _helpRequested(false), output_file_base("") {};
  virtual ~FriendsOrFoesApp();
  int get_num_clones(void) const { return num_clones; };
  bool get_is_rigid() const { return is_rigid; };
  const Clone &get_clone(int i) { return *clone_ptrs.at(i); };
  const std::string &get_output_file_base(void) { return output_file_base; }
  int get_maximum_time(void) { return maximum_time; }
  void printProperties(const std::string& base, std::ostream &ostrm) const;
  const Population *get_const_population_ptr(void) const { 
    return population_ptr; 
  }

protected:
	void initialize(Application& self);
	void uninitialize();
	void reinitialize(Application& self);

  virtual void defineOptions(OptionSet& options);

	void handleHelp(const std::string& name, const std::string& value);
	void handleConfig(const std::string& name, const std::string& value);
  void handleOutputFile(const std::string& name, const std::string& value);
  void handleNumClones(const std::string& name, const std::string& value);
  void handleRandomSeed(const std::string& name, const std::string& value);
  void parseVectorEntry(const std::string& name, const std::string& value, 
                        int &i, float &x);
  void parseMatrixEntry(const std::string& name, const std::string& value, 
                    int &i, int &j, float &x);
  void handleDefaultFitness(const std::string& name, const std::string& value);
  void handleSurvivalCoefficient(const std::string& name, const std::string& value);
  void handleReproductionCoefficient(const std::string& name, 
                                      const std::string& value);
  void handleCellwiseInteraction(const std::string& name, const std::string& value);
  void handleFitnesswiseInteraction(const std::string& name, const std::string& value);

	void displayHelp();
  virtual void set_values_from_config_with_defaults(void);
  virtual void set_num_clones(int new_num_clones);
  virtual HexPopulation *get_new_hex_population(void);
  virtual void create_population(void);
  virtual void destroy_population(void);
  void set_is_rigid(bool rigidity) { is_rigid = rigidity; };

  int num_clones;
  std::vector<Clone *> clone_ptrs;
  std::string output_file_base;
  bool is_rigid;
  int maximum_time;
  unsigned long random_seed;

	bool _helpRequested;
  Population *population_ptr;
};

#endif
