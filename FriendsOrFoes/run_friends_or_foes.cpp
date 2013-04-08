// run_friends_or_foes.cpp
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
#include "run_friends_or_foes.h"
#include "hex_replication_record.h"
#include <fstream>
#ifdef USING_MPI
#include <mpi.h>
#endif

RunFriendsOrFoesApp::~RunFriendsOrFoesApp() {
  if (cell_cycle_ptr != NULL) {
    delete cell_cycle_ptr;
  }
}

void RunFriendsOrFoesApp::create_cell_cycle(void) {
  if (is_rigid) {
    cell_cycle_ptr = new CellCycle(*population_ptr,
                                    make_hex_replication_record);
  }
}

int RunFriendsOrFoesApp::main(const vector<string>& args) {
  if (!_helpRequested)
  {
    set_values_from_config_with_defaults();
    string output_config = output_file_base + "_used.cnf";
    ofstream cnfstrm;
    cnfstrm.open(output_config.c_str(), ios::out);
    printProperties("fof", cnfstrm);
    cnfstrm.close();

    string output_file_name = output_file_base + (is_rigid ? ".hxg" : ".flx");
    ofstream ostrm;
    ostrm.open(output_file_name.c_str(), ios::out);
    try {
#ifdef USING_MPI
      // MPI-2 conformant MPI implementations are required to allow
      // applications to pass NULL for both the argc and argv arguments of
      // MPI_Init.  This is good since Poco::Application::init() has already
      // swallowed up the original argc and argv, and we would have to
      // reconstruct them....
      MPI_Init(NULL, NULL);
#endif
      Random::initialize(random_seed);
      create_population();
      create_cell_cycle();
      int generation;
      for (generation = 0; generation < maximum_time; ++generation) {
        if (population_ptr->get_num_cells() == 0) {
          break;
        }
        cell_cycle_ptr->run();
      }
    } catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        ostrm.close();
        Random::finalize();
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        ostrm.close();
        Random::finalize();
#ifdef USING_MPI
        MPI_Finalize();
#endif
    }
  }
  return Application::EXIT_OK;
}

POCO_APP_MAIN(FriendsOrFoesApp)
