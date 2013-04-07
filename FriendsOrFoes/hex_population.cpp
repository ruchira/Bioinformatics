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
#include "hex_population.h"
#include "hex_cell.h"
#include "probability.h"
#include <cassert>
#include <gsl/gsl_randist.h>

void HexPopulation::fill_field_with_clone(Clone &clone) {
  HexCell *hex_cell_ptr;
  int i;
  for (i = 0, hex_cell_ptr = hex_cell_grid; i < total_num_possible_cells; 
      ++i, ++hex_cell_ptr) {
    hex_cell_ptr->set_clone_ptr(&clone);
  }
}

void HexPopulation::make_focus_of_clone_at(Clone &clone, 
                                        int horiz_coord, int diag_coord) {
  HexCell &focus_cell = *cell_at(horiz_coord, diag_coord);
  focus_cell.set_clone_ptr(&clone);
}

void HexPopulation::make_focus_of_clone_in_middle(Clone &clone) {
  // Actually any horizontal coordinate is "the middle" since the grid is a
  // tube.  We'll use the middle coordinate anyway.
  make_focus_of_clone_at(clone, get_initial_width() / 2, 
                          get_initial_height() / 2);
}

void HexPopulation::make_focus_of_clone_at_random_spot(Clone &clone) {
  int horiz_coord = Random::rand_int_between_inclusive(0, 
                                                    get_initial_width() - 1);
  int diag_coord = Random::rand_int_between_inclusive(0, 
                                                    get_initial_height() - 1);
  
}

void HexPopulation::envivify(void) {
  HexCell *hex_cell_ptr;
  int i;
  for (i = 0, hex_cell_ptr = hex_cell_grid; i < total_num_possible_cells; 
      ++i, ++hex_cell_ptr) {
    hex_cell_ptr->envivify();
  }
}

void *HexPopulation::check_for_space_to_replicate(Cell &cell) {
  HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_dead_neighbors = fill_neighbor_buffer(*hex_cell_ptr, true);
  if (num_dead_neighbors > 0) {
    return neighbor_buffer[
                Random::rand_int_between_inclusive(1, num_dead_neighbors) - 1];
  }
  return NULL;
}

void HexPopulation::replicate(Cell &cell, void *space_specification,
                              ReplicationRecord &replication_record) {
  replication_record.record_as_mother(&cell);
  HexCell *daughter_cell_ptr = (HexCell *)space_specification;
  daughter_cell_ptr->set_clone_ptr(cell.get_clone_ptr());
  daughter_cell_ptr->envivify();
  replication_record.set_daughter0_ptr(&cell);
  replication_record.set_daughter1_ptr(daughter_cell_ptr);
}

/* There are six possible neighbors of each cell:
 5---4        ^ increasing row
 /\ /\       /
6--o--3     /
 \/ \/     /
 1---2     ----> increasing column
The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
making the grid topologically a tube like the esophagus
*/
bool HexPopulation::are_neighbors(const Cell &cell0, const Cell &cell1) const {
  const HexCell *hex_cell0_ptr = (HexCell *)(&cell0);
  const HexCell *hex_cell1_ptr = (HexCell *)(&cell1);
  int horiz_diff = (hex_cell1_ptr->get_horiz_coord() 
                    - hex_cell0_ptr->get_horiz_coord()) 
                    % get_initial_width();
  if (hex_cell1_ptr->get_diag_coord() == hex_cell0_ptr->get_diag_coord()) {
    if (horiz_diff == 1 || horiz_diff == -1) {
      return true;
    }
  } else if (hex_cell1_ptr->get_diag_coord() 
            == hex_cell0_ptr->get_diag_coord() - 1) {
    if (horiz_diff == 1 || horiz_diff == 0) {
      return true;
    }
  } else if (hex_cell1_ptr->get_diag_coord() 
            == hex_cell0_ptr->get_diag_coord() + 1) {
    if (horiz_diff == 0 || horiz_diff == -1) {
      return true;
    }
  }
  return false;
}

float HexPopulation::get_size_of_interface_between_cells(
                                    const Cell &cell0, const Cell &cell1) const {
  return (are_neighbors(cell0, cell1) ? 1.0 : 0.0);
}

HexPopulation::HexPopulation(int width, int height,
          function<void(Cell &, Cell &)> func) :
    Population(width, height, 1, func), 
    neighbor_permutation_size(6), 
    edge_neighbor_permutation_size(4) {
  assert(height >= 2);
  assert(width >= 3);
  total_num_possible_cells = width * height;
  // The grid consists of successive rows, each of which is contiguous in the
  // memory block.
  hex_cell_grid = new HexCell[total_num_possible_cells];
  global_permutation_size = total_num_possible_cells;
  neighbor_permutation_ptr = gsl_permutation_alloc(neighbor_permutation_size);
  edge_neighbor_permutation_ptr 
    = gsl_permutation_alloc(edge_neighbor_permutation_size);
  global_permutation_ptr = gsl_permutation_alloc(global_permutation_size);
  gsl_permutation_init(neighbor_permutation_ptr);
  gsl_permutation_init(edge_neighbor_permutation_ptr);
  gsl_permutation_init(global_permutation_ptr);
}

HexPopulation::~HexPopulation() {
  gsl_permutation_free(global_permutation_ptr);
  gsl_permutation_free(edge_neighbor_permutation_ptr);
  gsl_permutation_free(neighbor_permutation_ptr);
  delete hex_cell_grid;
}

HexCell *HexPopulation::cell_at(int horiz_coord, int diag_coord) {
  if (diag_coord >= 0 && diag_coord < get_initial_height()) {
    return &(hex_cell_grid[diag_coord * get_initial_width() 
                        + (horiz_coord % get_initial_width())]);
  }
  return NULL;
}

const HexCell *HexPopulation::const_cell_at(int horiz_coord, 
                                            int diag_coord) const {
  if (diag_coord >= 0 && diag_coord < get_initial_height()) {
    return &(hex_cell_grid[diag_coord * get_initial_width() 
                        + (horiz_coord % get_initial_width())]);
  }
  return NULL;
}

void HexPopulation::neighbor_reshuffle(void) {
  gsl_ran_shuffle(Random::get_gsl_rng_ptr(), neighbor_permutation_ptr->data,
                  neighbor_permutation_size, sizeof(size_t));
}

void HexPopulation::edge_neighbor_reshuffle(void) {
  gsl_ran_shuffle(Random::get_gsl_rng_ptr(), 
                  edge_neighbor_permutation_ptr->data,
                  edge_neighbor_permutation_size, sizeof(size_t));
}

void HexPopulation::global_reshuffle(void) const {
  gsl_ran_shuffle(Random::get_gsl_rng_ptr(), global_permutation_ptr->data,
                  global_permutation_size, sizeof(size_t));
}

void HexPopulation::map(function<void (Cell &)> proc, 
                  const Clone * filter_clone,
                  bool negate_filter) {
  int i;
  HexCell *hex_cell_ptr;
  if (filter_clone == NULL) {
    for (i = 0, hex_cell_ptr = hex_cell_grid; i < total_num_possible_cells; 
        ++i, ++hex_cell_ptr) {
      if (hex_cell_ptr->is_alive()) {
        proc(*hex_cell_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() != filter_clone) {
          proc(*hex_cell_ptr);
        }
      }
    } else {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() == filter_clone) {
          proc(*hex_cell_ptr);
        }
      }
    }
  }
} 

void HexPopulation::map_in_random_order(function<void (Cell &)> proc, 
                  const Clone * filter_clone,
                  bool negate_filter) {
  global_reshuffle();
  int i;
  if (filter_clone == NULL) {
    for (i = 0; i < total_num_possible_cells; ++i) {
      HexCell &hex_cell 
        = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
      if (hex_cell.is_alive()) {
        proc(hex_cell);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell = hex_cell_grid[
                              gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() && 
            hex_cell.get_clone_ptr() != filter_clone) {
          proc(hex_cell);
        }
      }
    } else {
      for (int i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell = hex_cell_grid[
                              gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() && 
            hex_cell.get_clone_ptr() == filter_clone) {
          proc(hex_cell);
        }
      }
    }
  }
} 

// This is like map except const.
void HexPopulation::const_map(function<void (const Cell &)> proc, 
                        const Clone * filter_clone,
                        bool negate_filter) const {
  int i;
  HexCell *hex_cell_ptr;
  if (filter_clone == NULL) {
    for (i = 0, hex_cell_ptr = hex_cell_grid; i < total_num_possible_cells; 
        ++i, ++hex_cell_ptr) {
      if (hex_cell_ptr->is_alive()) {
        proc(*hex_cell_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() != filter_clone) {
          proc(*hex_cell_ptr);
        }
      }
    } else {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() == filter_clone) {
          proc(*hex_cell_ptr);
        }
      }
    }
  }
}

// This is like map_in_random_order except const.
void HexPopulation::const_map_in_random_order(
                  function<void (const Cell &)> proc, 
                  const Clone * filter_clone,
                  bool negate_filter) const {
  global_reshuffle();
  int i;
  if (filter_clone == NULL) {
    for (i = 0; i < total_num_possible_cells; ++i) {
      HexCell &hex_cell 
        = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
      if (hex_cell.is_alive()) {
        proc(hex_cell);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell = hex_cell_grid[
                              gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() && 
            hex_cell.get_clone_ptr() != filter_clone) {
          proc(hex_cell);
        }
      }
    } else {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell = hex_cell_grid[
                              gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() && 
            hex_cell.get_clone_ptr() == filter_clone) {
          proc(hex_cell);
        }
      }
    }
  }
} 

void * HexPopulation::fold(function<void * (void *, Cell &)> func, 
                        void * initial_value_ptr, 
                        const Clone * filter_clone,
                        bool negate_filter) {
  void *data = initial_value_ptr;
  int i;
  HexCell *hex_cell_ptr;
  if (filter_clone == NULL) {
    for (i = 0, hex_cell_ptr = hex_cell_grid; i < total_num_possible_cells; 
        ++i, ++hex_cell_ptr) {
      if (hex_cell_ptr->is_alive()) {
        data = func(data, *hex_cell_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() != filter_clone) {
          data = func(data, *hex_cell_ptr);
        }
      }
    } else {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() == filter_clone) {
          data = func(data, *hex_cell_ptr);
        }
      }
    }
  }
  return data;
}

void * HexPopulation::fold_in_random_order(
                        function<void * (void *, Cell &)> func, 
                        void * initial_value_ptr, 
                        const Clone * filter_clone,
                        bool negate_filter) {
  global_reshuffle();
  void *data = initial_value_ptr;
  int i;
  if (filter_clone == NULL) {
    for (i = 0; i < total_num_possible_cells; ++i) {
      HexCell &hex_cell 
        = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
      if (hex_cell.is_alive()) {
        data = func(data, hex_cell);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell 
          = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() &&
            hex_cell.get_clone_ptr() != filter_clone) {
          data = func(data, hex_cell);
        }
      }
    } else {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell 
          = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() &&
            hex_cell.get_clone_ptr() == filter_clone) {
          data = func(data, hex_cell);
        }
      }
    }
  }
  return data;
}

void * HexPopulation::const_fold(function<void * (void *, const Cell &)> func, 
                              void * initial_value_ptr, 
                              const Clone * filter_clone,
                              bool negate_filter) const {
  void *data = initial_value_ptr;
  int i;
  HexCell *hex_cell_ptr;
  if (filter_clone == NULL) {
    for (i = 0, hex_cell_ptr = hex_cell_grid; i < total_num_possible_cells; 
        ++i, ++hex_cell_ptr) {
      if (hex_cell_ptr->is_alive()) {
        data = func(data, *hex_cell_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() != filter_clone) {
          data = func(data, *hex_cell_ptr);
        }
      }
    } else {
      for (i = 0, hex_cell_ptr = hex_cell_grid; 
          i < total_num_possible_cells; 
          ++i, ++hex_cell_ptr) {
        if (hex_cell_ptr->is_alive() && 
            hex_cell_ptr->get_clone_ptr() == filter_clone) {
          data = func(data, *hex_cell_ptr);
        }
      }
    }
  }
  return data;
}

void * HexPopulation::const_fold_in_random_order(
                        function<void * (void *, const Cell &)> func, 
                        void * initial_value_ptr, 
                        const Clone * filter_clone,
                        bool negate_filter) const {
  global_reshuffle();
  void *data = initial_value_ptr;
  int i;
  if (filter_clone == NULL) {
    for (i = 0; i < total_num_possible_cells; ++i) {
      HexCell &hex_cell 
        = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
      if (hex_cell.is_alive()) {
        data = func(data, hex_cell);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell 
          = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() &&
            hex_cell.get_clone_ptr() != filter_clone) {
          data = func(data, hex_cell);
        }
      }
    } else {
      for (i = 0; i < total_num_possible_cells; ++i) {
        HexCell &hex_cell 
          = hex_cell_grid[gsl_permutation_get(global_permutation_ptr, i)];
        if (hex_cell.is_alive() &&
            hex_cell.get_clone_ptr() == filter_clone) {
          data = func(data, hex_cell);
        }
      }
    }
  }
  return data;
}

void HexPopulation::map_neighbors(function<void (Cell &, Cell &)> proc,
                    Cell &cell,
                    const Clone * filter_neighbor_clone,
                    bool negate_filter) {
  HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_neighbors = fill_neighbor_buffer(*hex_cell_ptr);
  int i;
  HexCell **neighbor_ptrptr;
  if (filter_neighbor_clone == NULL) {
    for (i = 0, neighbor_ptrptr = neighbor_buffer; i < num_neighbors;
        ++i, ++neighbor_ptrptr) {
      if ((*neighbor_ptrptr)->is_alive()) {
        proc(cell, **neighbor_ptrptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0, neighbor_ptrptr = neighbor_buffer; i < num_neighbors;
          ++i, ++neighbor_ptrptr) {
        if ((*neighbor_ptrptr)->is_alive() && 
            (*neighbor_ptrptr)->get_clone_ptr() != filter_neighbor_clone) {
          proc(cell, **neighbor_ptrptr);
        }
      }
    } else {
      for (i = 0, neighbor_ptrptr = neighbor_buffer; i < num_neighbors;
          ++i, ++neighbor_ptrptr) {
        if ((*neighbor_ptrptr)->is_alive() && 
            (*neighbor_ptrptr)->get_clone_ptr() == filter_neighbor_clone) {
          proc(cell, **neighbor_ptrptr);
        }
      }
    }
  }
}

void HexPopulation::map_neighbors_in_random_order(
                    function<void (Cell &, Cell &)> proc,
                    Cell &cell,
                    const Clone * filter_neighbor_clone,
                    bool negate_filter) {
  HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_neighbors = fill_neighbor_buffer(*hex_cell_ptr);
  gsl_permutation *permutation_ptr;
  if (hex_cell_ptr->get_diag_coord() == 0 
      || hex_cell_ptr->get_diag_coord() == get_initial_height() - 1) {
    permutation_ptr = edge_neighbor_permutation_ptr;
    edge_neighbor_reshuffle();
  } else {
    permutation_ptr = neighbor_permutation_ptr;
    neighbor_reshuffle();
  }
  int i;
  if (filter_neighbor_clone == NULL) {
    for (i = 0; i < num_neighbors; ++i) {
      HexCell &neighbor 
        = *neighbor_buffer[gsl_permutation_get(permutation_ptr, i)];
      if (neighbor.is_alive()) {
        proc(cell, neighbor);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < num_neighbors; ++i) {
        HexCell &neighbor 
          = *neighbor_buffer[gsl_permutation_get(permutation_ptr, i)];
        if (neighbor.is_alive() && 
            neighbor.get_clone_ptr() != filter_neighbor_clone) {
          proc(cell, neighbor);
        }
      }
    } else {
      for (i = 0; i < num_neighbors; ++i) {
        HexCell &neighbor 
          = *neighbor_buffer[gsl_permutation_get(permutation_ptr, i)];
        if (neighbor.is_alive() && 
            neighbor.get_clone_ptr() == filter_neighbor_clone) {
          proc(cell, neighbor);
        }
      }
    }
  }
}

void HexPopulation::const_map_neighbors(
              function<void (const Cell &, const Cell &)> proc,
              const Cell &cell,
              const Clone * filter_neighbor_clone,
              bool negate_filter) const {
  const HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_neighbors = 6;
  if (hex_cell_ptr->get_diag_coord() == 0 
      || hex_cell_ptr->get_diag_coord() == get_initial_height() - 1) {
    num_neighbors = 4;
  }
  int i;
  if (filter_neighbor_clone == NULL) {
    for (i = 0; i < num_neighbors; ++i) {
      const HexCell *neighbor_ptr = const_get_neighbor(*hex_cell_ptr, i);
      if (neighbor_ptr->is_alive()) {
        proc(cell, *neighbor_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr = const_get_neighbor(*hex_cell_ptr, i);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() != filter_neighbor_clone) {
          proc(cell, *neighbor_ptr);
        }
      }
    } else {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr = const_get_neighbor(*hex_cell_ptr, i);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() == filter_neighbor_clone) {
          proc(cell, *neighbor_ptr);
        }
      }
    }
  }
}

void HexPopulation::const_map_neighbors_in_random_order(
              function<void (const Cell &, const Cell &)> proc,
              const Cell & cell,
              const Clone * filter_neighbor_clone,
              bool negate_filter) const {
  const HexCell *hex_cell_ptr = (const HexCell *)(&cell);
  int num_neighbors = 0;
  int permuted_neighbor_indices[6];
  if (hex_cell_ptr->get_diag_coord() == 0 
      || hex_cell_ptr->get_diag_coord() == get_initial_height() - 1) {
    num_neighbors = 4;
  } else {
    num_neighbors = 6;
  }
  int i, j;
  // Use the inside-out algorithm to get a an unbiased random permutation of
  // the neighbors.
  permuted_neighbor_indices[0] = 0;
  for (i = 1; i < num_neighbors; ++i) {
    j = Random::rand_int_between_inclusive(0, i);
    permuted_neighbor_indices[i] = permuted_neighbor_indices[j];
    permuted_neighbor_indices[j] = i;
  }
  if (filter_neighbor_clone == NULL) {
    for (i = 0; i < num_neighbors; ++i) {
      const HexCell *neighbor_ptr
        = const_get_neighbor(*hex_cell_ptr, permuted_neighbor_indices[i]);
      if (neighbor_ptr->is_alive()) {
        proc(cell, *neighbor_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr
          = const_get_neighbor(*hex_cell_ptr, permuted_neighbor_indices[i]);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() != filter_neighbor_clone) {
          proc(cell, *neighbor_ptr);
        }
      }
    } else {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr
          = const_get_neighbor(*hex_cell_ptr, permuted_neighbor_indices[i]);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() == filter_neighbor_clone) {
          proc(cell, *neighbor_ptr);
        }
      }
    }
  }
}

void * HexPopulation::fold_neighbors(
                    function<void * (void *, Cell &, Cell &)> func, 
                    Cell &cell,
                    void * initial_value_ptr, 
                    const Clone * filter_neighbor_clone,
                    bool negate_filter) {
  int i;
  HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_neighbors = fill_neighbor_buffer(*hex_cell_ptr);
  void *data = initial_value_ptr;
  HexCell **neighbor_ptrptr;
  if (filter_neighbor_clone == NULL) {
    for (i = 0, neighbor_ptrptr = neighbor_buffer; 
        i < num_neighbors;
        ++i, ++neighbor_ptrptr) {
      if ((*neighbor_ptrptr)->is_alive()) {
        data = func(data, cell, **neighbor_ptrptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0, neighbor_ptrptr = neighbor_buffer; 
          i < num_neighbors;
          ++i, ++neighbor_ptrptr) {
        if ((*neighbor_ptrptr)->is_alive() && 
            (*neighbor_ptrptr)->get_clone_ptr() != filter_neighbor_clone) {
          data = func(data, cell, **neighbor_ptrptr);
        }
      }
    } else {
      for (i = 0, neighbor_ptrptr = neighbor_buffer; 
          i < num_neighbors;
          ++i, ++neighbor_ptrptr) {
        if ((*neighbor_ptrptr)->is_alive() && 
            (*neighbor_ptrptr)->get_clone_ptr() == filter_neighbor_clone) {
          data = func(data, cell, **neighbor_ptrptr);
        }
      }
    }
  }
  return data;
}

void * HexPopulation::fold_neighbors_in_random_order(
                    function<void * (void *, Cell &, Cell &)> func, 
                    Cell &cell,
                    void * initial_value_ptr, 
                    const Clone * filter_neighbor_clone,
                    bool negate_filter) {
  HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_neighbors = fill_neighbor_buffer(*hex_cell_ptr);
  gsl_permutation *permutation_ptr;
  if (hex_cell_ptr->get_diag_coord() == 0 
      || hex_cell_ptr->get_diag_coord() == get_initial_height() - 1) {
    permutation_ptr = edge_neighbor_permutation_ptr;
    edge_neighbor_reshuffle();
  } else {
    permutation_ptr = neighbor_permutation_ptr;
    neighbor_reshuffle();
  }
  void *data = initial_value_ptr;
  int i;
  if (filter_neighbor_clone == NULL) {
    for (i = 0; i < num_neighbors; ++i) {
      HexCell &neighbor
        = *neighbor_buffer[gsl_permutation_get(permutation_ptr, i)];
      if (neighbor.is_alive()) {
        data = func(data, cell, neighbor);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < num_neighbors; ++i) {
        HexCell &neighbor
          = *neighbor_buffer[gsl_permutation_get(permutation_ptr, i)];
        if (neighbor.is_alive() && 
            neighbor.get_clone_ptr() != filter_neighbor_clone) {
          data = func(data, cell, neighbor);
        }
      }
    } else {
      for (i = 0; i < num_neighbors; ++i) {
        HexCell &neighbor
          = *neighbor_buffer[gsl_permutation_get(permutation_ptr, i)];
        if (neighbor.is_alive() && 
            neighbor.get_clone_ptr() == filter_neighbor_clone) {
          data = func(data, cell, neighbor);
        }
      }
    }
  }
  return data;
}

void * HexPopulation::const_fold_neighbors(
        function<void * (void *, const Cell &, const Cell &)> func, 
        const Cell &cell,
        void * initial_value_ptr, 
        const Clone * filter_neighbor_clone,
        bool negate_filter) const {
  const HexCell *hex_cell_ptr = (HexCell *)(&cell);
  int num_neighbors = 6;
  if (hex_cell_ptr->get_diag_coord() == 0 
      || hex_cell_ptr->get_diag_coord() == get_initial_height() - 1) {
    num_neighbors = 4;
  }
  void *data = initial_value_ptr;
  int i;
  if (filter_neighbor_clone == NULL) {
    for (i = 0; i < num_neighbors; ++i) {
      const HexCell *neighbor_ptr = const_get_neighbor(*hex_cell_ptr, i);
      if (neighbor_ptr->is_alive()) {
        data = func(data, cell, *neighbor_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr = const_get_neighbor(*hex_cell_ptr, i);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() != filter_neighbor_clone) {
          data = func(data, cell, *neighbor_ptr);
        }
      }
    } else {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr = const_get_neighbor(*hex_cell_ptr, i);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() == filter_neighbor_clone) {
          data = func(data, cell, *neighbor_ptr);
        }
      }
    }
  }
  return data;
}

void * HexPopulation::const_fold_neighbors_in_random_order(
        function<void * (void *, const Cell &, const Cell &)> func, 
        const Cell &cell,
        void * initial_value_ptr, 
        const Clone * filter_neighbor_clone,
        bool negate_filter) const {
  void *data = initial_value_ptr;
  const HexCell *hex_cell_ptr = (const HexCell *)(&cell);
  int num_neighbors = 0;
  int permuted_neighbor_indices[6];
  if (hex_cell_ptr->get_diag_coord() == 0 
      || hex_cell_ptr->get_diag_coord() == get_initial_height() - 1) {
    num_neighbors = 4;
  } else {
    num_neighbors = 6;
  }
  int i, j;
  // Use the inside-out algorithm to get a an unbiased random permutation of
  // the neighbors.
  permuted_neighbor_indices[0] = 0;
  for (i = 1; i < num_neighbors; ++i) {
    j = Random::rand_int_between_inclusive(0, i);
    permuted_neighbor_indices[i] = permuted_neighbor_indices[j];
    permuted_neighbor_indices[j] = i;
  }
  if (filter_neighbor_clone == NULL) {
    for (i = 0; i < num_neighbors; ++i) {
      const HexCell *neighbor_ptr
        = const_get_neighbor(*hex_cell_ptr, permuted_neighbor_indices[i]);
      if (neighbor_ptr->is_alive()) {
        data = func(data, cell, *neighbor_ptr);
      }
    }
  } else {
    if (negate_filter) {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr
          = const_get_neighbor(*hex_cell_ptr, permuted_neighbor_indices[i]);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() != filter_neighbor_clone) {
          data = func(data, cell, *neighbor_ptr);
        }
      }
    } else {
      for (i = 0; i < num_neighbors; ++i) {
        const HexCell *neighbor_ptr
          = const_get_neighbor(*hex_cell_ptr, permuted_neighbor_indices[i]);
        if (neighbor_ptr->is_alive() && 
            neighbor_ptr->get_clone_ptr() == filter_neighbor_clone) {
          data = func(data, cell, *neighbor_ptr);
        }
      }
    }
  }
  return data;
}

void HexPopulation::clear_neighbor_buffer() {
  for (int i = 0; i < 6; ++i) {
    neighbor_buffer[i] = NULL;
  }
}

/* There are six possible neighbors of each cell:
 5---4        ^ increasing row
 /\ /\       /
6--o--3     /
 \/ \/     /
 1---2     ----> increasing column
The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
making the grid topologically a tube like the esophagus
*/
HexCell *HexPopulation::get_neighbor(HexCell &cell, int i) {
  HexCell *cell_ptr = &cell;
  switch(i) {
    case 1:
      if (cell.get_diag_coord() > 0) {
        // Go down to the previous row (i.e., back a whole row of cells in the
        // memory block).
        cell_ptr -= get_initial_width();
        return cell_ptr;
      } else {
        // This cell is at the bottom edge and has no neighbor 1
        return NULL;
      }
      break;
    case 2:
      if (cell.get_diag_coord() > 0) {
        // Go down to the previous row (i.e., back a whole row of cells in the
        // memory block).
        cell_ptr -= get_initial_width();
        // Go right by one (diagonal) column.
        if (cell.get_horiz_coord() == get_initial_width() - 1) {
          // Since we are at (diagonal) column (get_initial_width() - 1), to go
          // right by one column is the same as to go left by
          // (get_initial_width() - 1) columns.
          cell_ptr -= (get_initial_width() - 1);
        } else {
          ++cell_ptr; 
        }
        return cell_ptr;
      } else {
        // This cell is at the bottom edge and has no neighbor 1
        return NULL;
      }
      break;
    case 3:
      // Go right by one (diagonal) column.
      if (cell.get_horiz_coord() == get_initial_width() - 1) {
        // Since we are at (diagonal) column (get_initial_width() - 1), to go
        // right by one column is the same as to go left by
        // (get_initial_width() - 1) columns.
        cell_ptr -= (get_initial_width() - 1);
      } else {
        ++cell_ptr;
      }
      return cell_ptr;
      break;
    case 4:
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // Go up to the next row (i.e. forward a whole row of cells in the
        // memory block.
        cell_ptr += get_initial_width();
        return cell_ptr;
      } else {
        // This cell is at the top edge and has no neighbor 4
        return NULL;
      }
      break;
    case 5:
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // Go up to the next row (i.e. forward a whole row of cells in the
        // memory block.
        cell_ptr += get_initial_width();
        // Go left by one (diagonal) column.
        if (cell.get_horiz_coord() == 0) {
          // Since we are at (diagonal) column 0, to go left by one column is
          // the same as to go right by (get_initial_width() - 1) columns.
          cell_ptr += (get_initial_width() - 1);
        } else {
          --cell_ptr;
        }
        return cell_ptr;
      } else {
        // This cell is at the top edge and has no neighbor 4
        return NULL;
      }
      break;
    case 6:
      // Go left by one (diagonal) column.
      if (cell.get_horiz_coord() == 0) {
        // Since we are at (diagonal) column 0, to go left by one column is
        // the same as to go right by (get_initial_width() - 1) columns.
        cell_ptr += (get_initial_width() - 1);
      } else {
        --cell_ptr;
      }
      return cell_ptr;
      break;
    default:
      return NULL;
  };
}

/* There are six possible neighbors of each cell:
 5---4        ^ increasing row
 /\ /\       /
6--o--3     /
 \/ \/     /
 1---2     ----> increasing column
The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
making the grid topologically a tube like the esophagus
*/
const HexCell *HexPopulation::const_get_neighbor(const HexCell &cell, 
                                                  int i) const {
  const HexCell *cell_ptr = &cell;
  switch(i) {
    case 1:
      if (cell.get_diag_coord() > 0) {
        // Go down to the previous row (i.e., back a whole row of cells in the
        // memory block).
        cell_ptr -= get_initial_width();
        return cell_ptr;
      } else {
        // This cell is at the bottom edge and has no neighbor 1
        return NULL;
      }
      break;
    case 2:
      if (cell.get_diag_coord() > 0) {
        // Go down to the previous row (i.e., back a whole row of cells in the
        // memory block).
        cell_ptr -= get_initial_width();
        // Go right by one (diagonal) column.
        if (cell.get_horiz_coord() == get_initial_width() - 1) {
          // Since we are at (diagonal) column (get_initial_width() - 1), to go
          // right by one column is the same as to go left by
          // (get_initial_width() - 1) columns.
          cell_ptr -= (get_initial_width() - 1);
        } else {
          ++cell_ptr; 
        }
        return cell_ptr;
      } else {
        // This cell is at the bottom edge and has no neighbor 1
        return NULL;
      }
      break;
    case 3:
      // Go right by one (diagonal) column.
      if (cell.get_horiz_coord() == get_initial_width() - 1) {
        // Since we are at (diagonal) column (get_initial_width() - 1), to go
        // right by one column is the same as to go left by
        // (get_initial_width() - 1) columns.
        cell_ptr -= (get_initial_width() - 1);
      } else {
        ++cell_ptr;
      }
      return cell_ptr;
      break;
    case 4:
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // Go up to the next row (i.e. forward a whole row of cells in the
        // memory block.
        cell_ptr += get_initial_width();
        return cell_ptr;
      } else {
        // This cell is at the top edge and has no neighbor 4
        return NULL;
      }
      break;
    case 5:
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // Go up to the next row (i.e. forward a whole row of cells in the
        // memory block.
        cell_ptr += get_initial_width();
        // Go left by one (diagonal) column.
        if (cell.get_horiz_coord() == 0) {
          // Since we are at (diagonal) column 0, to go left by one column is
          // the same as to go right by (get_initial_width() - 1) columns.
          cell_ptr += (get_initial_width() - 1);
        } else {
          --cell_ptr;
        }
        return cell_ptr;
      } else {
        // This cell is at the top edge and has no neighbor 4
        return NULL;
      }
      break;
    case 6:
      // Go left by one (diagonal) column.
      if (cell.get_horiz_coord() == 0) {
        // Since we are at (diagonal) column 0, to go left by one column is
        // the same as to go right by (get_initial_width() - 1) columns.
        cell_ptr += (get_initial_width() - 1);
      } else {
        --cell_ptr;
      }
      return cell_ptr;
      break;
    default:
      return NULL;
  };
}

/* There are six possible neighbors of each cell:
 5---4        ^ increasing row
 /\ /\       /
6--o--3     /
 \/ \/     /
 1---2     ----> increasing column
The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
making the grid topologically a tube like the esophagus
*/
int HexPopulation::fill_neighbor_buffer(HexCell &cell, bool require_dead) {
  clear_neighbor_buffer();
  int num_neighbors = 0;
  HexCell **neighbor_ptrptr;
  neighbor_ptrptr = neighbor_buffer;
  HexCell *cell_ptr = &cell;
  if (require_dead) {
    if (cell.get_horiz_coord() == 0) {
      // We are at the leftmost column
      if (cell.get_diag_coord() > 0) {
        // neighbor 1: go down to the previous row
        cell_ptr -= get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 2: go right by one (diagonal) column
        ++cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 3: go back up to the original row
        cell_ptr += get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      } else {
        // This cell is at the bottom, so has no neighbors 1 and 2
        // neighbor 3: go right by one (diagonal) column
        ++cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      }
      // now we are at neighbor 3
      // go back left by one column, ending up at the original cell
      --cell_ptr;
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // neighbor 4: go up to the next row
        cell_ptr += get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 5: go left by one (diagonal) column
        // Since we are at column 0, this is the same as going right by
        // (get_initial_width() - 1) (diagonal) columns.
        cell_ptr += (get_initial_width() - 1);
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 6: go back down to the original row
        cell_ptr -= get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      } else {
        // This cell is at the top, so has no neighbors 4 and 5
        // neighbor 6: go left by one (diagonal) column
        // Since we are at column 0, this is the same as going right by
        // (get_initial_width() - 1) (diagonal) columns.
        cell_ptr += (get_initial_width() - 1);
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      }
    } else if (cell.get_horiz_coord() == get_initial_width() - 1) {
      // We are at the rightmost column
      if (cell.get_diag_coord() > 0) {
        // neighbor 1: go down to the previous row
        cell_ptr -= get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 2: go right by one (diagonal) column
        // Since we are at (diagonal) column (get_initial_width() - 1), to go right
        // by one column is the same as to go left by (get_initial_width() - 1)
        // columns.
        cell_ptr -= (get_initial_width() - 1);
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 3: go back up to the original row
        cell_ptr += get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      } else {
        // This cell is at the bottom, so has no neighbors 1 and 2
        // neighbor 3: go right by one (diagonal) column
        // Since we are at (diagonal) column (get_initial_width() - 1), to go right
        // by one column is the same as to go left by (get_initial_width() - 1)
        // columns.
        cell_ptr -= (get_initial_width() - 1);
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      }
      // now we are at neighbor 3, which is in the leftmost column
      // go back left by one column, ending up at the original cell
      // Since we are at (diagonal) column 0, to go left by one column is the
      // same as to go right by (get_initial_width() - 1) columns.
      cell_ptr += (get_initial_width() - 1);
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // neighbor 4: go up to the next row
        cell_ptr += get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 5: go left by one (diagonal) column
        --cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 6: go back down to the original row
        cell_ptr -= get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      } else {
        // This cell is at the top, so has no neighbors 4 and 5
        // neighbor 6: go left by one (diagonal) column
        --cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      }
    } else {
      // We are neither at the leftmost column nor the rightmost column
      if (cell.get_diag_coord() > 0) {
        // neighbor 1: go down to the previous row
        cell_ptr -= get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 2: go right by one (diagonal) column
        ++cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 3: go back up to the original row
        cell_ptr += get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      } else {
        // This cell is at the bottom, so has no neighbors 1 and 2
        // neighbor 3: go right by one (diagonal) column
        ++cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      }
      // go back left by one column, ending up at the original cell
      --cell_ptr;
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // neighbor 4: go up to the next row
        cell_ptr += get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 5: go left by one (diagonal) column
        --cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
        // neighbor 6: go back down to the original row
        cell_ptr -= get_initial_width();
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      } else {
        // This cell is at the top, so has no neighbors 4 and 5
        // neighbor 6: go left by one (diagonal) column
        --cell_ptr;
        if (!cell_ptr->is_alive()) {
          *neighbor_ptrptr++ = cell_ptr;
          ++num_neighbors;
        }
      }
    }
  } else {
    // require_dead is false
    if (cell.get_horiz_coord() == 0) {
      // We are at the leftmost column
      if (cell.get_diag_coord() > 0) {
        // neighbor 1: go down to the previous row
        cell_ptr -= get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 2: go right by one (diagonal) column
        ++cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 3: go back up to the original row
        cell_ptr += get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      } else {
        // This cell is at the bottom, so has no neighbors 1 and 2
        // neighbor 3: go right by one (diagonal) column
        ++cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      }
      // now we are at neighbor 3
      // go back left by one column, ending up at the original cell
      --cell_ptr;
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // neighbor 4: go up to the next row
        cell_ptr += get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 5: go left by one (diagonal) column
        // Since we are at column 0, this is the same as going right by
        // (get_initial_width() - 1) (diagonal) columns.
        cell_ptr += (get_initial_width() - 1);
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 6: go back down to the original row
        cell_ptr -= get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      } else {
        // This cell is at the top, so has no neighbors 4 and 5
        // neighbor 6: go left by one (diagonal) column
        // Since we are at column 0, this is the same as going right by
        // (get_initial_width() - 1) (diagonal) columns.
        cell_ptr += (get_initial_width() - 1);
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      }
    } else if (cell.get_horiz_coord() == get_initial_width() - 1) {
      // We are at the rightmost column
      if (cell.get_diag_coord() > 0) {
        // neighbor 1: go down to the previous row
        cell_ptr -= get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 2: go right by one (diagonal) column
        // Since we are at (diagonal) column (get_initial_width() - 1), to go right
        // by one column is the same as to go left by (get_initial_width() - 1)
        // columns.
        cell_ptr -= (get_initial_width() - 1);
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 3: go back up to the original row
        cell_ptr += get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      } else {
        // This cell is at the bottom, so has no neighbors 1 and 2
        // neighbor 3: go right by one (diagonal) column
        // Since we are at (diagonal) column (get_initial_width() - 1), to go right
        // by one column is the same as to go left by (get_initial_width() - 1)
        // columns.
        cell_ptr -= (get_initial_width() - 1);
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      }
      // now we are at neighbor 3, which is in the leftmost column
      // go back left by one column, ending up at the original cell
      // Since we are at (diagonal) column 0, to go left by one column is the
      // same as to go right by (get_initial_width() - 1) columns.
      cell_ptr += (get_initial_width() - 1);
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // neighbor 4: go up to the next row
        cell_ptr += get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 5: go left by one (diagonal) column
        --cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 6: go back down to the original row
        cell_ptr -= get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      } else {
        // This cell is at the top, so has no neighbors 4 and 5
        // neighbor 6: go left by one (diagonal) column
        --cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      }
    } else {
      // We are neither at the leftmost column nor the rightmost column
      if (cell.get_diag_coord() > 0) {
        // neighbor 1: go down to the previous row
        cell_ptr -= get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 2: go right by one (diagonal) column
        ++cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 3: go back up to the original row
        cell_ptr += get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      } else {
        // This cell is at the bottom, so has no neighbors 1 and 2
        // neighbor 3: go right by one (diagonal) column
        ++cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      }
      // go back left by one column, ending up at the original cell
      --cell_ptr;
      if (cell.get_diag_coord() < get_initial_height() - 1) {
        // neighbor 4: go up to the next row
        cell_ptr += get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 5: go left by one (diagonal) column
        --cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
        // neighbor 6: go back down to the original row
        cell_ptr -= get_initial_width();
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      } else {
        // This cell is at the top, so has no neighbors 4 and 5
        // neighbor 6: go left by one (diagonal) column
        --cell_ptr;
        *neighbor_ptrptr++ = cell_ptr;
        ++num_neighbors;
      }
    }
  }
  return num_neighbors;
}
