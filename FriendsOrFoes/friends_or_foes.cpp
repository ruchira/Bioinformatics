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

#include <iostream>
#include <limits>
#include <math.h>
#include <stdlib.h>

#include "Poco/Util/Application.h"
#include "Poco/Util/Option.h"
#include "Poco/Util/OptionSet.h"
#include "Poco/Util/IntValidator.h"
#include "Poco/Util/HelpFormatter.h"
#include "Poco/Util/AbstractConfiguration.h"
#include "Poco/AutoPtr.h"

#include "clone.h"
#include "friends_or_foes.h"

using Poco::Util::Application;
using Poco::Util::Option;
using Poco::Util::OptionSet;
using Poco::Util::HelpFormatter;
using Poco::Util::AbstractConfiguration;
using Poco::Util::OptionCallback;
using Poco::Util::IntValidator;
using Poco::AutoPtr;

const float square_root_of_three = sqrt(3.0);
const pair<float,float> unit_diagonal = make_pair(0.5,
                                                  0.5 * square_root_of_three);

class FriendsOrFoesApp: public Application {
public:
  FriendsOrFoesApp(): _helpRequested(false) {
  }

protected:
	void initialize(Application& self) {
    num_clones = 0;  
		loadConfiguration(); // load default configuration files, if present
		Application::initialize(self);
	}

	void uninitialize() {
		Application::uninitialize();
	}

	void reinitialize(Application& self) {
		Application::reinitialize(self);
	}

	void defineOptions(OptionSet& options) {
		Application::defineOptions(options);

		options.addOption(
			Option("help", "?", "display help information on command line arguments")
				.required(false)
				.repeatable(false)
				.callback(OptionCallback<FriendsOrFoesApp>(this, 
                                              &FriendsOrFoesApp::handleHelp)));

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
      Option("maximum_time", "T", "maximum time step of the simulation")
        .required(false)
        .repeatable(false)
        .argument("value", true)
        .validator(new IntValidator(1,numeric_limits<int>::max()))
        .binding("fof.maximum_time"));

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
      Option("num_clones", "n", "number of clones")
        .required(false)
        .repeatable(false)
        .argument("value", true)
        .validator(new IntValidator(1,numeric_limits<int>::max()))
        .callback(OptionCallback<FriendsOrFoesApp>(this,
                                          &FriendsOrFoesApp::handleNumClones)));

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

		options.addOption(
			Option("config-file", "f", "load configuration data from a file")
				.required(false)
				.repeatable(true)
				.argument("file")
				.callback(OptionCallback<FriendsOrFoesApp>(this, &FriendsOrFoesApp::handleConfig)));

	}
	
	void handleHelp(const string& name, const string& value) {
		_helpRequested = true;
		displayHelp();
		stopOptionsProcessing();
	}
	
  void handleNumClones(const string& name, const string& value) {
    config().setString("fof." + name, value);
    set_num_clones(strtol(value.c_str(), NULL, 10));
  }

  void parseVectorEntry(const string& name, const string& value, 
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

  void parseMatrixEntry(const string& name, const string& value, 
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
        col.assign(value, col_pos + 1, entry_pos - col_pos);
        i = strtol(row.c_str(), NULL, 10);
        j = strtol(col.c_str(), NULL, 10);
        x = strtof(entry.c_str(), NULL);
        config().setString(("fof." + name + '_' + row + '_' + col), entry);
      }
    }
  }

  void handleDefaultFitness(const string& name, const string& value) {
    int i;
    float x;
    parseVectorEntry(name, value, i, x);
    clone_ptrs.at(i)->set_default_fitness(x);
  }

  void handleSurvivalCoefficient(const string& name, 
                                  const string& value) {
    int i;
    float x;
    parseVectorEntry(name, value, i, x);
    clone_ptrs.at(i)->set_survival_coefficient(x);
  }

  void handleReproductionCoefficient(const string& name, 
                                      const string& value) {
    int i;
    float x;
    parseVectorEntry(name, value, i, x);
    clone_ptrs.at(i)->set_reproduction_coefficient(x);
  }

  void handleCellwiseInteraction(const string& name, 
                                  const string& value) {
    
    int i, j;
    float x;
    parseMatrixEntry(name, value, i, j, x);
    clone_ptrs.at(i)->set_constant_effect_of_clone(clone_ptrs.at(j), x);
  }

  void handleFitnesswiseInteraction(const string& name, 
                                    const string& value) {
    int i, j;
    float x;
    parseMatrixEntry(name, value, i, j, x);
    clone_ptrs.at(i)->set_linear_effect_of_clone(clone_ptrs.at(j), x);
  }

	void handleConfig(const string& name, const string& value) {
		loadConfiguration(value);
	}
		
	void displayHelp() {
		HelpFormatter helpFormatter(options());
		helpFormatter.setCommand(commandName());
		helpFormatter.setUsage("OPTIONS");
		helpFormatter.setHeader("Application to simulate coexistence, competition and cooperation in an epithelium.\nNB: Clones are numbered from 0.");
    helpFormatter.setWidth(256);
		helpFormatter.format(cout);
	}
	
	void defineProperty(const string& def) {
		string name;
		string value;
		string::size_type pos = def.find('=');
		if (pos != string::npos)
		{
			name.assign(def, 0, pos);
			value.assign(def, pos + 1, def.length() - pos);
		}
		else name = def;
		config().setString(name, value);
	}

	int main(const vector<string>& args) {
		if (!_helpRequested)
		{
			logger().information("Arguments to main():");
			for (vector<string>::const_iterator it = args.begin(); it != args.end(); ++it)
			{
				logger().information(*it);
			}
			logger().information("Application properties:");
			printProperties("");
      try {
        cout << "Under development!" << endl;
      } catch(exception& e) {
          cerr << "error: " << e.what() << "\n";
          return 1;
      }
      catch(...) {
          cerr << "Exception of unknown type!\n";
      }
		}
		return Application::EXIT_OK;
	}
	
	void printProperties(const string& base) {
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
			}
		}
		else
		{
			for (AbstractConfiguration::Keys::const_iterator it = keys.begin(); it != keys.end(); ++it)
			{
				string fullKey = base;
				if (!fullKey.empty()) fullKey += '.';
				fullKey.append(*it);
				printProperties(fullKey);
			}
		}
	}

  int get_num_clones() {
    return num_clones;
  }

private:
	bool _helpRequested;
  int num_clones;
  vector<Clone *> clone_ptrs;

  void set_num_clones(int new_num_clones) {
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
};

POCO_APP_MAIN(FriendsOrFoesApp)
