// hex_population.h
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
#ifndef HEX_POPULATION_H
#define HEX_POPULATION_H
#include "population.h"
#include "hex_cell.h"
#include <gsl/gsl_permutation.h>

class HexPopulation : public Population {
  public:
    HexPopulation(int width, int height,
          function<void(Cell &, Cell &)> func
            =neighbor_affect_cell_through_affine_function);
    virtual ~HexPopulation();
    // This fills the grid with cells of the specified clone.
    void fill_field_with_clone(Clone &clone);
    // This puts a cell of the specified clone at the specified spot.
    void make_focus_of_clone_at(Clone &clone, int horiz_coord, int diag_coord);
    void make_focus_of_clone_in_middle(Clone &clone);
    void make_focus_of_clone_at_random_spot(Clone &clone);
    // This envivifies all cells.
    void envivify(void);

    /*
     * Methods to act globally or locally on the population of cells
     */

    virtual void map(function<void (Cell &)> proc, 
                      const Clone * filter_clone = NULL,
                      bool negate_filter = false);
    virtual void map_in_random_order(function<void (Cell &)> proc,
                                      const Clone * filter_clone = NULL,
                                      bool negate_filter = false);
    virtual void const_map(function<void (const Cell &)> proc, 
                            const Clone * filter_clone = NULL,
                            bool negate_filter = false) const;
    virtual void const_map_in_random_order(function<void (const Cell &)> proc,
                                    const Clone * filter_clone = NULL,
                                    bool negate_filter = false) const;
    virtual void * fold(function<void * (void *, Cell &)> func, 
                        void * initial_value_ptr, 
                        const Clone * filter_clone = NULL,
                        bool negate_filter = false);
    virtual void * fold_in_random_order(function<void * (void *, Cell &)> func,
                                        void * initial_value_ptr,
                                        const Clone * filter_clone = NULL,
                                        bool negate_filter = false);
    virtual void * const_fold(function<void * (void *, const Cell &)> func, 
                              void * initial_value_ptr, 
                              const Clone * filter_clone = NULL,
                              bool negate_filter = false) const;
    virtual void * const_fold_in_random_order(
                                  function<void * (void *, const Cell &)> func,
                                  void * initial_value_ptr,
                                  const Clone * filter_clone = NULL,
                                  bool negate_filter = false) const;
    virtual void map_neighbors(function<void (Cell &, Cell &)> proc,
                                Cell &cell,
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false);
    virtual void map_neighbors_in_random_order(
                                function<void (Cell &, Cell &)> proc,
                                Cell &cell,
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false);
    virtual void const_map_neighbors(
                          function<void (const Cell &, const Cell &)> proc,
                          const Cell &cell,
                          const Clone * filter_neighbor_clone = NULL,
                          bool negate_filter = false) const;
    virtual void const_map_neighbors_in_random_order(
                          function<void (const Cell &, const Cell &)> proc,
                          const Cell & cell_ptr,
                          const Clone * filter_neighbor_clone = NULL,
                          bool negate_filter = false) const;
    virtual void * fold_neighbors(
                                function<void * (void *, Cell &, Cell &)> func, 
                                Cell &cell,
                                void * initial_value_ptr, 
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false);
    virtual void * fold_neighbors_in_random_order(
                                function<void * (void *, Cell &, Cell &)> func, 
                                Cell &cell,
                                void * initial_value_ptr, 
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false);
    virtual void * const_fold_neighbors(
                    function<void * (void *, const Cell &, const Cell &)> func, 
                    const Cell &cell,
                    void * initial_value_ptr, 
                    const Clone * filter_neighbor_clone = NULL,
                    bool negate_filter = false) const;
    virtual void * const_fold_neighbors_in_random_order(
                    function<void * (void *, const Cell &, const Cell &)> func, 
                    const Cell &cell,
                    void * initial_value_ptr, 
                    const Clone * filter_neighbor_clone = NULL,
                    bool negate_filter = false) const;
    /*
     * End of methods to act globally or locally on the population of cells
     */

    HexCell * cell_at(int horiz_coord, int diag_coord);
    const HexCell * const_cell_at(int horiz_coord, int diag_coord) const;

    virtual bool are_neighbors(const Cell &cell0, const Cell &cell1) const;

    /* There are six possible neighbors of each cell:
     5---4        ^ increasing row
     /\ /\       /
    6--o--3     /
     \/ \/     /
     1---2     ----> increasing column
    The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
    making the grid topologically a tube like the esophagus
    */
    // The following two functions each return a reference to the ith neighbor
    // of the cell, using the numbering above, or NULL if it doesn't exist
    // (either because i is not a number from 1 to 6, or because the cell is at
    // the top or bottom edge of the grid and doesn't have a neighbor in that
    // direction).
    HexCell* get_neighbor(HexCell &cell, int i);
    const HexCell* const_get_neighbor(const HexCell &cell, int i) const;

    // This returns 0.0 if the cells are not neighbors.
    virtual float get_size_of_interface_between_cells(
                                                  const Cell &cell0, 
                                                  const Cell &cell1) const;
    virtual float get_median_distance_from_clone_to_clone(const Clone &clone,
                                    const Clone &neighbor_clone) const;
    // This returns NULL if there is no space to replicate, or a description of
    // the available space if there is.
    virtual void *check_for_space_to_replicate(Cell &cell);
    // This replicates cell into the specified space and records the
    // replication in replication_record.
    virtual void replicate(Cell &cell, void *space_specification,
                          ReplicationRecord &replication_record);

  private:
    HexCell *hex_cell_grid;
    HexCell *neighbor_buffer[6];
    void clear_neighbor_buffer();
    int fill_neighbor_buffer(HexCell &cell, bool require_dead = false);
    const size_t neighbor_permutation_size;
    const size_t edge_neighbor_permutation_size;
    size_t global_permutation_size;
    gsl_permutation *neighbor_permutation_ptr;
    gsl_permutation *edge_neighbor_permutation_ptr;
    gsl_permutation *global_permutation_ptr;
    void neighbor_reshuffle(void);
    void edge_neighbor_reshuffle(void);
    void global_reshuffle(void) const;
    int total_num_possible_cells;
};

#endif
