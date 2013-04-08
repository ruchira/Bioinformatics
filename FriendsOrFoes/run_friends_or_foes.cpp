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
#include <iostream>
#include <fstream>
#ifdef USING_MPI
#include <mpi.h>
#endif

using namespace friends_or_foes;

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

struct DataForWritingReplicationRecord {
  DataForWritingReplicationRecord(HexPopulationEvent &event, ostream &strm) :
    hex_population_event(event), ostrm(strm) {};
  HexPopulationEvent &hex_population_event;
  ostream &ostrm;
};

void *write_replication_to_stream(void *data, const ReplicationRecord &record)
{
  DataForWritingReplicationRecord *pStruct 
    = (DataForWritingReplicationRecord *)data;
  const HexReplicationRecord *hex_record_ptr 
        = (HexReplicationRecord *)(&record);
  HexCellProto *hex_cell_proto0_ptr 
    = pStruct->hex_population_event.mutable_cell0();
  HexCellProto *hex_cell_proto1_ptr
    = pStruct->hex_population_event.mutable_cell1();
  hex_cell_proto0_ptr->set_horiz_coord(
        hex_record_ptr->get_hex_mother_ptr()->get_horiz_coord());
  hex_cell_proto0_ptr->set_diag_coord(
        hex_record_ptr->get_hex_mother_ptr()->get_diag_coord());
  hex_cell_proto1_ptr->set_horiz_coord(
        hex_record_ptr->get_hex_daughter1_ptr()->get_horiz_coord());
  hex_cell_proto1_ptr->set_diag_coord(
        hex_record_ptr->get_hex_daughter1_ptr()->get_diag_coord());
  pStruct->hex_population_event.SerializeToOstream(&pStruct->ostrm);
  return data;
}

void RunFriendsOrFoesApp::write_hex_cell_cycle_run(ostream &ostrm) {
  HexCellProto *hex_cell_proto0_ptr;
  HexPopulationEvent hex_population_event;
  int i, j;
  hex_population_event.Clear();
  hex_population_event.set_type(HexPopulationEvent::kill);
  hex_cell_proto0_ptr = hex_population_event.mutable_cell0();
  for (i = 0; i < num_clones; ++i) {
    const vector<const Cell *> *killed_cells 
      = cell_cycle_ptr->get_killed_cells_of_clone(*clone_ptrs[i]);
    for (j = 0; j < killed_cells->size(); ++j) {
      hex_cell_proto0_ptr->set_horiz_coord(
        ((const HexCell *)killed_cells->at(j))->get_horiz_coord());
      hex_cell_proto0_ptr->set_diag_coord(
        ((const HexCell *)killed_cells->at(j))->get_diag_coord());
      hex_population_event.SerializeToOstream(&ostrm);
    }
  }
  hex_population_event.set_type(HexPopulationEvent::replicate);
  DataForWritingReplicationRecord data(hex_population_event, ostrm);
  cell_cycle_ptr->const_fold_replication_records(write_replication_to_stream, 
                                                  &data);
  hex_population_event.clear_cell0();
  hex_population_event.clear_cell1();
  hex_population_event.set_type(HexPopulationEvent::stop);
  hex_population_event.SerializeToOstream(&ostrm);
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
        write_hex_cell_cycle_run(ostrm);
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
