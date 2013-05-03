// see_friends_or_foes.h
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
#ifndef SEE_FRIENDS_OR_FOES_H
#define SEE_FRIENDS_OR_FOES_H
#include "friends_or_foes.h"
#include "hex_friends_or_foes_history.pb.h"
#include "hex_cell.h"
#include "hex_replication_record.h"

class SeeFriendsOrFoesApp : public FriendsOrFoesApp {
  public:
    SeeFriendsOrFoesApp() : FriendsOrFoesApp() {};
    virtual ~SeeFriendsOrFoesApp() {};
    virtual void defineOptions(OptionSet& options);
    int main(const std::vector<std::string>& args);
		// This will read HexPopulationEvents from the stream until an event with
		// the type STOP is encountered.  It returns whether or not it read any
		// events (including the STOP event).
		bool read_and_replay_hex_cell_cycle_run(std::istream &istrm);
	protected:
		virtual HexPopulation *get_new_hex_population(void);
		virtual void visualize(int generation);
    virtual void set_values_from_config_with_defaults(void);
  private:
    float max_fitness_ever;
    friends_or_foes::HexPopulationEvent hex_population_event;
		HexReplicationRecord hex_replication_record;
    HexCell* get_hex_cell_ptr_of_hex_cell_proto(
                                            const friends_or_foes::HexCellProto&);
};

#endif
