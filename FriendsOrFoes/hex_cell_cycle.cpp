// hex_cell_cycle.cpp
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
#include "hex_cell_cycle.h"

void HexCellCycle::run(void) {
	live_neighbors_of_dead_cells.clear();
	CellCycle::run();
}

void *add_live_neighbor_of_dead_cell(void *data, Cell &cell, Cell &neighbor) {
	HexCellCycle *hex_cell_cycle_ptr = (HexCellCycle *)data;
	if (neighbor.is_alive()) {
		hex_cell_cycle_ptr->live_neighbors_of_dead_cells_set.insert(&neighbor);
	}
	return data;
}

void HexCellCycle::try_survival(Cell &cell) {
	CellCycle::try_survival(cell);
	if (!cell.is_alive()) {
		live_neighbors_of_dead_cells_set.erase(&cell);
		population.fold_neighbors(add_live_neighbor_of_dead_cell, cell, this);
	}
}

void HexCellCycle::global_try_reproduction(void) {
	live_neighbors_of_dead_cells.clear();
	std::set<Cell *>::const_iterator iter;
	for (iter = live_neighbors_of_dead_cells_set.begin(); 
			iter != live_neighbors_of_dead_cells_set.end(); ++iter) {
		live_neighbors_of_dead_cells.push_back(*iter);
	}
	// The cells in the set went into the vector in an arbitrary order.  However,
	// this order is not guaranteed to be a uniform random ordering.
	// Use the Fisher-Yates algorithm to get a uniform random ordering.
	int i, j;
	Cell *temp;
	for (i = live_neighbors_of_dead_cells.size() - 1; i >= 1; --i) {
		j = Random::rand_int_between_inclusive(0, i);
		temp = live_neighbors_of_dead_cells[j];
		live_neighbors_of_dead_cells[j] = live_neighbors_of_dead_cells[i];
		live_neighbors_of_dead_cells[i] = temp;
	}
	// Now check the reproduction of these cells, in the permuted order.
	// No cell can reproduce unless there is space next to it (i.e., a dead
	// cell), so these are all the cells we need to check.
	for (i = 0; i < live_neighbors_of_dead_cells.size(); ++i) {
		try_reproduction(*live_neighbors_of_dead_cells[i]);
	}
}
