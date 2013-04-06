// friends_or_foes.cpp
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
#include "friends_or_foes.h"

#include <limits>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "Poco/Util/IntValidator.h"
#include "Poco/Util/HelpFormatter.h"
#include "Poco/Util/AbstractConfiguration.h"
#include "Poco/AutoPtr.h"

#include "friends_or_foes.h"
#include "replication_record.h"

using Poco::Util::Application;
using Poco::Util::Option;
using Poco::Util::OptionSet;
using Poco::Util::HelpFormatter;
using Poco::Util::AbstractConfiguration;
using Poco::Util::OptionCallback;
using Poco::Util::IntValidator;
using Poco::AutoPtr;

void FriendsOrFoesApp::initialize(Application& self) {
  set_num_clones(2);
  loadConfiguration(); // load default configuration files, if present
  Application::initialize(self);
}

void FriendsOrFoesApp::uninitialize() {
  Application::uninitialize();
}

void FriendsOrFoesApp::reinitialize(Application& self) {
  Application::reinitialize(self);
}

void FriendsOrFoesApp::defineOptions(OptionSet& options) {
  Application::defineOptions(options);

  options.addOption(
    Option("help", "?", "display help information on command line arguments")
      .required(false)
      .repeatable(false)
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                                            &FriendsOrFoesApp::handleHelp)));

  options.addOption(
    Option("output_file_base", "o", "base of the output file")
      .required(false)
      .repeatable(false)
      .argument("file")
      .callback(OptionCallback<FriendsOrFoesApp>(this,
                                      &FriendsOrFoesApp::handleOutputFile)));

  options.addOption(
    Option("config_file", "c", "load configuration data from a file")
      .required(false)
      .repeatable(true)
      .argument("file")
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                                          &FriendsOrFoesApp::handleConfig)));

  options.addOption(
    Option("num_clones", "n", "number of clones")
      .required(false)
      .repeatable(false)
      .argument("value", true)
      .validator(new IntValidator(1,numeric_limits<int>::max()))
      .callback(OptionCallback<FriendsOrFoesApp>(this,
                                        &FriendsOrFoesApp::handleNumClones)));

  options.addOption(
    Option("rigid", "r", "use a rigid grid")
      .required(false)
      .repeatable(false)
      .group("grid_type")
      .binding("fof.rigid"));

  options.addOption(
    Option("Voronoi", "V", "use a Voronoi grid")
      .required(false)
      .repeatable(false)
      .group("grid_type")
      .binding("fof.voronoi"));

  options.addOption(
    Option("maximum_time", "T", "maximum time step of the simulation")
      .required(false)
      .repeatable(false)
      .argument("value", true)
      .validator(new IntValidator(1,numeric_limits<int>::max()))
      .binding("fof.maximum_time"));

  options.addOption(
    Option("random_seed", "r", "seed for the random number generator")
      .required(false)
      .repeatable(false)
      .argument("value", true)
      .validator(new IntValidator(1,numeric_limits<unsigned long>::max()))
      .binding("fof.random_seed"));

  options.addOption(
    Option("width", "w", "width in cells of the initial grid")
      .required(false)
      .repeatable(false)
      .argument("value", true)
      .validator(new IntValidator(1,numeric_limits<int>::max()))
      .binding("fof.width"));

  options.addOption(
    Option("height", "h", "height in cells of the initial grid")
      .required(false)
      .repeatable(false)
      .argument("value", true)
      .validator(new IntValidator(1,numeric_limits<int>::max()))
      .binding("fof.height"));

  options.addOption(
    Option("minimum_divisible_area", "a", 
            "minimum area at which a Voronoi cell can divide")
    .required(false)
    .repeatable(false)
    .argument("value", true)
    .binding("fof.minimum_divisible_area"));

  options.addOption(
    Option("default_fitness", "D", "default fitness of a clone")
      .required(false)
      .repeatable(true)
      .argument("clone_number:fitness", true)
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                                &FriendsOrFoesApp::handleDefaultFitness)));

  options.addOption(
    Option("survival_coefficient", "S", 
        "fitness multiplier giving survival probability of a cell of a clone")
      .required(false)
      .repeatable(true)
      .argument("clone_number:coefficient", true)
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                              &FriendsOrFoesApp::handleSurvivalCoefficient)));

  options.addOption(
    Option("reproduction_coefficient", "R", 
    "fitness multiplier giving reproduction probability of a cell of a clone")
      .required(false)
      .repeatable(true)
      .argument("clone_number:coefficient", true)
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                          &FriendsOrFoesApp::handleReproductionCoefficient)));

  options.addOption(
    Option("cellwise_interaction", "E", 
    "amount to add to cell of clone 0's fitness per each neighboring cell of clone 1")
    .required(false)
    .repeatable(true)
    .argument("clone_number_0,clone_number_1:increment", true)
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                              &FriendsOrFoesApp::handleCellwiseInteraction)));

  options.addOption(
    Option("fitnesswise_interaction", "F", 
    "fraction of the fitness of each neighboring cell of clone 1 to add to a cell of clone 0's fitness")
    .required(false)
    .repeatable(true)
    .argument("clone_number_0,clone_number_1:fraction", true)
      .callback(OptionCallback<FriendsOrFoesApp>(this, 
                          &FriendsOrFoesApp::handleFitnesswiseInteraction)));

}
	
void FriendsOrFoesApp::handleHelp(const string& name, const string& value) {
  _helpRequested = true;
  displayHelp();
  stopOptionsProcessing();
}

void FriendsOrFoesApp::handleConfig(const string& name, const string& value) {
  loadConfiguration(value, PRIO_APPLICATION);
}

void FriendsOrFoesApp::handleOutputFile(const string& name, const string& value) 
{
  output_file_name = value;
  config().setString("fof.output_file_base", value);
}

void FriendsOrFoesApp::handleNumClones(const string& name, const string& value) {
  config().setString("fof.num_clones", value);
  set_num_clones(strtol(value.c_str(), NULL, 10));
}

void FriendsOrFoesApp::parseVectorEntry(const string& name, const string& value, 
                                        int &i, float &x) {
  string index;
  string entry;
  string::size_type pos = value.find(':');
  if (pos != string::npos) {
    index.assign(value, 0, pos);
    entry.assign(value, pos + 1, value.length() - pos);
    i = strtol(index.c_str(), NULL, 10);
    x = strtof(entry.c_str(), NULL);
    config().setString(("fof." + name + '_' + index), entry);
  }
}

void FriendsOrFoesApp::parseMatrixEntry(const string& name, const string& value, 
                                        int &i, int &j, float &x) {
  string row;
  string col;
  string entry;
  string::size_type entry_pos = value.find(':');
  if (entry_pos != string::npos) {
    entry.assign(value, entry_pos + 1, value.length() - entry_pos);
    string::size_type col_pos = value.find(',');
    if (col_pos != string::npos) {
      row.assign(value, 0, col_pos);
      col.assign(value, col_pos + 1, entry_pos - col_pos - 1);
      i = strtol(row.c_str(), NULL, 10);
      j = strtol(col.c_str(), NULL, 10);
      x = strtof(entry.c_str(), NULL);
      config().setString(("fof." + name + '_' + row + '_' + col), entry);
    }
  }
}

void FriendsOrFoesApp::handleDefaultFitness(const string& name, 
                                            const string& value) {
  int i;
  float x;
  parseVectorEntry(name, value, i, x);
  clone_ptrs.at(i)->set_default_fitness(x);
} 

void FriendsOrFoesApp::handleSurvivalCoefficient(const string& name, 
                                                  const string& value) {
  int i;
  float x;
  parseVectorEntry(name, value, i, x);
  clone_ptrs.at(i)->set_survival_coefficient(x);
}

void FriendsOrFoesApp::handleReproductionCoefficient(const string& name, 
                                                      const string& value) {
  int i;
  float x;
  parseVectorEntry(name, value, i, x);
  clone_ptrs.at(i)->set_reproduction_coefficient(x);
}

void FriendsOrFoesApp::handleCellwiseInteraction(const string& name, 
                                                  const string& value) {
  
  int i, j;
  float x;
  parseMatrixEntry(name, value, i, j, x);
  clone_ptrs.at(i)->set_constant_effect_of_clone(clone_ptrs.at(j), x);
}

void FriendsOrFoesApp::handleFitnesswiseInteraction(const string& name, 
                                                    const string& value) {
  int i, j;
  float x;
  parseMatrixEntry(name, value, i, j, x);
  clone_ptrs.at(i)->set_linear_effect_of_clone(clone_ptrs.at(j), x);
}
  
void FriendsOrFoesApp::displayHelp() {
  HelpFormatter helpFormatter(options());
  helpFormatter.setCommand(commandName());
  helpFormatter.setUsage("OPTIONS");
  helpFormatter.setHeader("Application to simulate coexistence, competition and cooperation in an epithelium.\nNB: Clones are numbered from 0.");
  helpFormatter.setWidth(256);
  helpFormatter.format(cout);
}
	
void FriendsOrFoesApp::set_values_from_config_with_defaults(void) {
  if (!config().hasProperty("fof.num_clones")) {
    stringstream strvalue;
    strvalue << num_clones;
    config().setString("fof.num_clones", strvalue.str());
  }
  if (!config().hasProperty("fof.random_seed")) {
    config().setString("fof.random_seed", "12345");
  }
  random_seed = strtoul(config().getString("fof.random_seed").c_str(), NULL, 
                        10);
  if (!config().hasProperty("fof.rigid") 
      && !config().hasProperty("fof.voronoi")) {
    config().setString("fof.rigid", "");
  }
  if (config().hasProperty("fof.voronoi")) {
    is_rigid = false;
  } else {
    is_rigid = true;
  }
  if (!config().hasProperty("fof.output_file_base")) {
    config().setString("fof.output_file_base", "friends_or_foes_output");
  }
  if (!config().hasProperty("fof.maximum_time")) {
    config().setString("fof.maximum_time", "5000");
  }
  maximum_time = strtol(config().getString("fof.maximum_time").c_str(), NULL, 
                        10);
  output_file_name = config().getString("fof.output_file_base") 
                    + (is_rigid ? ".hxg" : ".flx");
  if (!config().hasProperty("fof.width")) {
    config().setString("fof.width", "1500");
  }
  if (!config().hasProperty("fof.height")) {
    config().setString("fof.height", "2000");
  }
  if (!is_rigid && !config().hasProperty("fof.minimum_divisible_area")) {
    config().setString("fof.minimum_divisible_area", "5.0");
  }
  for (int i = 0; i < num_clones; ++i) {
    stringstream default_fitness_name;
    default_fitness_name << "fof.default_fitness_" << i;
    if (!config().hasProperty(default_fitness_name.str())) {
      stringstream strvalue;
      strvalue << clone_ptrs[i]->get_default_fitness();
      config().setString(default_fitness_name.str(), strvalue.str());
    }
    stringstream survival_coefficient_name;
    survival_coefficient_name << "fof.survival_coefficient_" << i;
    if (!config().hasProperty(survival_coefficient_name.str())) {
      stringstream strvalue;
      strvalue << clone_ptrs[i]->get_survival_coefficient();
      config().setString(survival_coefficient_name.str(), strvalue.str());
    }
    stringstream reproduction_coefficient_name;
    reproduction_coefficient_name << "fof.reproduction_coefficient_" << i;
    if (!config().hasProperty(reproduction_coefficient_name.str())) {
      stringstream strvalue;
      strvalue << clone_ptrs[i]->get_reproduction_coefficient();
      config().setString(reproduction_coefficient_name.str(), strvalue.str());
    }
    for (int j = 0; j < num_clones; ++j) {
      stringstream cellwise_interaction_name;
      cellwise_interaction_name << "fof.cellwise_interaction_"  << i 
                                                                << "_" << j;
      if (!config().hasProperty(cellwise_interaction_name.str())) {
        stringstream strvalue;
        strvalue<< clone_ptrs[i]->get_constant_effect_of_clone(clone_ptrs[j]);
        config().setString(cellwise_interaction_name.str(), strvalue.str());
      }
      stringstream fitnesswise_interaction_name;
      fitnesswise_interaction_name << "fof.fitnesswise_interaction_"  << i 
                                                                << "_" << j;
      if (!config().hasProperty(fitnesswise_interaction_name.str())) {
        stringstream strvalue;
        strvalue << clone_ptrs[i]->get_linear_effect_of_clone(clone_ptrs[j]);
        config().setString(fitnesswise_interaction_name.str(),
                          strvalue.str());
      }
    }
  }
}
	
void FriendsOrFoesApp::printProperties(const string& base, ostream &ostrm) const 
{
  AbstractConfiguration::Keys keys;
  config().keys(base, keys);
  if (keys.empty())
  {
    if (config().hasProperty(base))
    {
      string msg;
      msg.append(base);
      msg.append(" = ");
      msg.append(config().getString(base));
      logger().information(msg);
      ostrm << msg << endl;
    }
  }
  else
  {
    for (AbstractConfiguration::Keys::const_iterator it = keys.begin(); 
        it != keys.end(); ++it)
    {
      string fullKey = base;
      if (!fullKey.empty()) fullKey += '.';
      fullKey.append(*it);
      printProperties(fullKey, ostrm);
    }
  }
}

void FriendsOrFoesApp::set_num_clones(int new_num_clones) {
  if (new_num_clones != num_clones) {
    Clone * clone_ptr;
    while (clone_ptrs.size() > 0) {
      clone_ptr = clone_ptrs.back();
      delete clone_ptr;
      clone_ptrs.pop_back();
    }
    num_clones = new_num_clones;
    for (int i = 0; i < num_clones; ++i) {
      clone_ptr = new Clone();
      clone_ptrs.push_back(clone_ptr);
    }
  }
}
