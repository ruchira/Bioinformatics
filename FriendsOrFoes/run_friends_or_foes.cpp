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
#include "hex_cell_cycle.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#ifdef USING_MPI
#include <mpi.h>
#endif
#include <google/protobuf/message_lite.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace friends_or_foes;
using namespace google::protobuf;
using namespace google::protobuf::io;

RunFriendsOrFoesApp::~RunFriendsOrFoesApp() {
  if (cell_cycle_ptr != NULL) {
    delete cell_cycle_ptr;
  }
}

void RunFriendsOrFoesApp::create_cell_cycle(void) {
  if (is_rigid) {
    cell_cycle_ptr = new HexCellCycle(*population_ptr);
  }
}

struct DataForWritingReplicationRecord {
  DataForWritingReplicationRecord(HexPopulationEvent &event, 
                                  CodedOutputStream &output) :
      hex_population_event(event), coded_output(output) {
    first = true;
  };
  HexPopulationEvent &hex_population_event;
  CodedOutputStream &coded_output;
  bool first;
  uint32 byte_size;
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
  if (pStruct->first) {
    // All replication events have the same size, so we only compute the size
    // the first time through the loop.
    pStruct->byte_size = pStruct->hex_population_event.ByteSize();
    pStruct->first = false;
  }
  pStruct->coded_output.WriteVarint32(pStruct->byte_size);
  pStruct->hex_population_event.SerializeWithCachedSizes(&pStruct->coded_output);
  return data;
}

void RunFriendsOrFoesApp::write_hex_cell_cycle_run(
                      CodedOutputStream &coded_output) {
  HexCellProto *hex_cell_proto0_ptr;
  HexPopulationEvent hex_population_event;
  int i, j;
  hex_population_event.Clear();
  hex_population_event.set_type(HexPopulationEvent::kill);
  hex_cell_proto0_ptr = hex_population_event.mutable_cell0();
  bool first = true;
  uint32 byte_size;
  for (i = 0; i < num_clones; ++i) {
    const std::vector<const Cell *> *killed_cells 
      = cell_cycle_ptr->get_killed_cells_of_clone(*clone_ptrs[i]);
    if (killed_cells != NULL) {
      for (j = 0; j < killed_cells->size(); ++j) {
        hex_cell_proto0_ptr->set_horiz_coord(
          ((const HexCell *)killed_cells->at(j))->get_horiz_coord());
        hex_cell_proto0_ptr->set_diag_coord(
          ((const HexCell *)killed_cells->at(j))->get_diag_coord());
        // All kill events have the same size, so we only compute the size the
        // first time through the loop.
        if (first) {
          byte_size = hex_population_event.ByteSize();
          first = false;
        }
        coded_output.WriteVarint32(byte_size);
        hex_population_event.SerializeWithCachedSizes(&coded_output);
      }
    }
  }
  hex_population_event.set_type(HexPopulationEvent::replicate);
  DataForWritingReplicationRecord data(hex_population_event, coded_output);
  cell_cycle_ptr->const_fold_replication_records(write_replication_to_stream, 
                                                  &data);
  hex_population_event.clear_cell0();
  hex_population_event.clear_cell1();
  hex_population_event.set_type(HexPopulationEvent::stop);
  byte_size = hex_population_event.ByteSize();
  coded_output.WriteVarint32(byte_size);
  hex_population_event.SerializeWithCachedSizes(&coded_output);
}

int RunFriendsOrFoesApp::main(const std::vector<std::string>& args) {
  if (!_helpRequested)
  {
    set_values_from_config_with_defaults();
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    std::string output_config = output_file_base + ".properties";
    std::ofstream cnfstrm;
    cnfstrm.open(output_config.c_str(), std::ios::out);
    printProperties("fof", cnfstrm);
    cnfstrm.close();

    std::string output_file_name 
      = output_file_base + (is_rigid ? ".hxg" : ".flx");
    int fd = open(output_file_name.c_str(), O_WRONLY);
    ZeroCopyOutputStream *raw_output = new FileOutputStream(fd);
    CodedOutputStream* coded_output = new CodedOutputStream(raw_output);
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
      int magic_number = 1769;
      coded_output->WriteLittleEndian32(magic_number);
      create_population();
      create_cell_cycle();
      int generation;
      for (generation = 0; generation < maximum_time; ++generation) {
        if (population_ptr->get_num_cells() == 0) {
          break;
        }
        cell_cycle_ptr->run();
        write_hex_cell_cycle_run(*coded_output);
      }
      std::string max_fitness_file_name = output_file_base + ".max";
      std::ofstream maxstrm;
      maxstrm.open(max_fitness_file_name.c_str(), std::ios::out);
      maxstrm << population_ptr->get_max_fitness_ever();
      maxstrm.close();
      destroy_population();
      delete coded_output;
      delete raw_output;
      close(fd);
      Random::finalize();
#ifdef USING_MPI
      MPI_Finalize();
#endif
      google::protobuf::ShutdownProtobufLibrary();
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        destroy_population();
        delete coded_output;
        delete raw_output;
        close(fd);
        Random::finalize();
#ifdef USING_MPI
        MPI_Finalize();
#endif
        google::protobuf::ShutdownProtobufLibrary();
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
        destroy_population();
        delete coded_output;
        delete raw_output;
        close(fd);
        google::protobuf::ShutdownProtobufLibrary();
        Random::finalize();
#ifdef USING_MPI
        MPI_Finalize();
#endif
    }
  }
  return Application::EXIT_OK;
}

POCO_APP_MAIN(RunFriendsOrFoesApp)
