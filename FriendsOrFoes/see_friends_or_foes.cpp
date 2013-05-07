// see_friends_or_foes.cpp
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
#include "see_friends_or_foes.h"
#include "visualize_hex_population.h"
#include "hexagon_renderings.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include "Poco/Util/IntValidator.h"


using namespace friends_or_foes;

void SeeFriendsOrFoesApp::defineOptions(OptionSet& options) {
  FriendsOrFoesApp::defineOptions(options);
  options.addOption(
    Option("magnification", "m", "magnification step at which to render grid (0-4)")
      .required(false)
      .repeatable(false)
      .argument("value", true)
      .validator(new IntValidator(0,4))
      .binding("fof.magnification"));
}

HexPopulation *SeeFriendsOrFoesApp::get_new_hex_population(void) {
  VisualizeHexPopulation *hex_population_ptr = new VisualizeHexPopulation(
      strtol(config().getString("fof.width").c_str(), NULL, 10),
      strtol(config().getString("fof.height").c_str(), NULL, 10),
      clone_ptrs,
      get_the_hexagon_rendering(
        strtol(config().getString("fof.magnification").c_str(), NULL, 10),
        num_clones),
      max_fitness_ever);
  return hex_population_ptr;
}

void SeeFriendsOrFoesApp::visualize(int generation) {
	std::ostringstream os;
	os << output_file_base << '_' 
      << std::setw(4) << std::setfill('0') << generation
			<< ".bmp";
	((VisualizeHexPopulation *)population_ptr)->visualize(os.str().c_str());
}

HexCell *SeeFriendsOrFoesApp::get_hex_cell_ptr_of_hex_cell_proto(
                                              const HexCellProto &hex_cell_proto) {
	return ((HexPopulation *)population_ptr)->cell_at(hex_cell_proto.horiz_coord(),
                                                    hex_cell_proto.diag_coord());
}

bool SeeFriendsOrFoesApp::read_and_replay_hex_cell_cycle_run(std::istream &istrm) {
	bool result = false;
	bool read_event = false;
	do {
		hex_population_event.Clear();
		read_event = hex_population_event.ParseFromIstream(&istrm);
		if (read_event) {
			result = true;
			if (hex_population_event.type() == HexPopulationEvent::kill) {
				population_ptr->kill(*get_hex_cell_ptr_of_hex_cell_proto(
																							hex_population_event.cell0()));
			}
			if (hex_population_event.type() == HexPopulationEvent::replicate) {
				population_ptr->replicate(*get_hex_cell_ptr_of_hex_cell_proto(
																							hex_population_event.cell0()),
																	get_hex_cell_ptr_of_hex_cell_proto(
																							hex_population_event.cell1()),
																	hex_replication_record);
			}
		}
	} while (read_event && 
					hex_population_event.type() != HexPopulationEvent::stop);
	return result;
}

void SeeFriendsOrFoesApp::set_values_from_config_with_defaults(void) {
  FriendsOrFoesApp::set_values_from_config_with_defaults();
  if (!config().hasProperty("fof.magnification")) {
    config().setString("fof.magnification", "0");
  }
}

int SeeFriendsOrFoesApp::main(const std::vector<std::string>& args) {
	if (!_helpRequested)
	{
    set_values_from_config_with_defaults();
    std::string max_fitness_file_name = output_file_base + ".max";
    std::ifstream max_fitness_strm;
    max_fitness_strm.open(max_fitness_file_name.c_str(), std::ios::in);
    max_fitness_strm >> max_fitness_ever;
    max_fitness_strm.close();
    std::string input_file_name = output_file_base + (is_rigid ? ".hxg" : ".flx");
    std::ifstream istrm;
    istrm.open(input_file_name.c_str(), std::ios::in);
    try {
      allegro_init();
      set_color_depth(8);
      create_population();
      readkey();
      int generation;
      for (generation = 0; generation < maximum_time; ++generation) {
				visualize(generation);
        readkey();
        if (!read_and_replay_hex_cell_cycle_run(istrm)) {
					break;
				}
      }
      visualize(generation);
      readkey();
			istrm.close();
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        istrm.close();
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
        istrm.close();
    }
  }
  return Application::EXIT_OK;
}

POCO_APP_MAIN(SeeFriendsOrFoesApp)
