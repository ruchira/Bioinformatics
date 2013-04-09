// hex_cell_cycle.h
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
#ifndef HEX_CELL_CYCLE_H
#define HEX_CELL_CYCLE_H
#include "cell_cycle.h"
#include "hex_replication_record.h"
#include <set>

class HexCellCycle : public CellCycle {
	public:
		HexCellCycle(Population &population) :
			CellCycle(population, make_hex_replication_record) {};
		virtual ~HexCellCycle() {};
		virtual void run(void);
	protected:
		virtual void try_survival(Cell &cell);
		virtual void global_try_reproduction(void);
	private:
		std::set<Cell *> live_neighbors_of_dead_cells_set;
		// We keep these here so the capacity will grow to the needed size and
		// persist between runs (even though the vector itself will be cleared).
		vector<Cell *> live_neighbors_of_dead_cells;
		friend void *add_live_neighbor_of_dead_cell(void *data, Cell &cell,
																								Cell &neighbor);
};
#endif
