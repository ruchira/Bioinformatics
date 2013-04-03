// population.cpp
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
#include "population.h"

void * count_cells(void * counter_ptr, const Cell &cell) {
  (*( (int *) counter_ptr ) )++;
  return counter_ptr;
}

int Population::get_num_cells(Clone *filter_clone) const {
  int result = 0;
  const_fold(&count_cells, &result, filter_clone);
  return result;
}

void * sum_volume_of_cells(void * cumulative_volume_ptr, const Cell &cell) {
  (*( (float *) cumulative_volume_ptr ) ) += cell.get_volume();
  return cumulative_volume_ptr;
}

float Population::get_volume(Clone *filter_clone) const {
  float result = 0;
  const_fold(&sum_volume_of_cells, &result, filter_clone);
  return result;
}

void * has_neighbors(void * flag_ptr, const Cell &cell, const Cell &neighbor) {
  *( (bool *) flag_ptr ) = true;
  return flag_ptr;
}

bool Population::is_adjacent_to_clone(const Cell &cell, 
                                              const Clone &clone) const {
  bool result = false;
  const_fold_neighbors(&has_neighbors, cell, &result, &clone);
}

struct AdjacentCountResult {
  AdjacentCountResult(const Population &a_population, const Clone &a_clone) : 
    population(a_population), neighbor_clone(a_clone), 
    num_cells_adjacent_to_neighbor_clone(0) {};
  const Population &population;
  const Clone &neighbor_clone;
  int num_cells_adjacent_to_neighbor_clone;
};

void * count_cells_adjacent_to_clone(void *count_ptr, const Cell &cell) {
  AdjacentCountResult *result_ptr = (AdjacentCountResult *)count_ptr;
  if (result_ptr->population.is_adjacent_to_clone(cell,
                                                result_ptr->neighbor_clone)) {
    ++(result_ptr->num_cells_adjacent_to_neighbor_clone);
  }
  return count_ptr;
}

int Population::get_num_cells_of_clone_adjacent_to_clone(
                                            const Clone &clone, 
                                            const Clone &neighbor_clone) const {
  AdjacentCountResult result = AdjacentCountResult(*this, neighbor_clone);
  const_fold(&count_cells_adjacent_to_clone, &result, &clone);
  return result.num_cells_adjacent_to_neighbor_clone;
}

struct BoundarySizeResult {
  BoundarySizeResult(const Population &a_population, const Clone &a_clone) :
      population(a_population), clone(a_clone), boundary_size(0.0) {};
  const Population &population;
  const Clone &clone;
  float boundary_size;
};

void *sum_cell_neighbor_interface(void *data, const Cell &cell0, 
                                              const Cell&cell1) {
  BoundarySizeResult *result_ptr = (BoundarySizeResult *)data;
  result_ptr->boundary_size += 
      result_ptr->population.get_size_of_interface_between_cells(cell0, cell1);
  return result_ptr;
}

void *sum_interface_of_cell_to_other_clones(void *data, const Cell &cell) {
  BoundarySizeResult *result_ptr = (BoundarySizeResult *)data;
  result_ptr->population.const_fold_neighbors(
    &sum_cell_neighbor_interface, cell, result_ptr, &(result_ptr->clone),
    true);
  return data;
}

float Population::get_size_of_boundary_of_clone(
                                                  const Clone &clone) const {
  
  BoundarySizeResult result = BoundarySizeResult(*this, clone);
  const_fold(&sum_interface_of_cell_to_other_clones, &result, &clone);
  return result.boundary_size;
}

void *update_fitness(void *data, Cell &cell) {
  if (cell.is_alive()) {
    cell.reset_fitness();
    Population *population_ptr;
    population_ptr = (Population *)data;
    population_ptr->map_neighbors(
                        population_ptr->get_neighbor_affect_cell_func(), cell);
  }
  return data;
}

void Population::update_all_fitnesses(void) {
  fold(update_fitness, this);
}

