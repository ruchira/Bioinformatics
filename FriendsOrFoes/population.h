// population.h
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
#ifndef POPULATION_H
#define POPULATION_H
#include <tr1/functional>
#include <map>
#include <vector>
#include <set>
#include "replication_record.h"

using std::tr1::function;

extern void neighbor_affect_cell_through_affine_function(Cell &cell, 
                                                          Cell &neighbor);

class Population {
  public:
    Population(int width, int height, int depth=1,
          function<void(Cell &, Cell &)> func
            =neighbor_affect_cell_through_affine_function) : 
        initial_width(width), initial_height(height), initial_depth(depth),
        neighbor_affect_cell_func(func), max_fitness_ever(0.0) {};
    virtual ~Population();
    int get_initial_width(void) const { return initial_width; };
    int get_initial_height(void) const { return initial_height; };
    int get_initial_depth(void) const { return initial_depth; };

    /*
     * Methods to act globally or locally on the population of cells
     */

    // For every cell in the population belonging to the clone filter_clone (or
    // every cell, if filter_clone is NULL; or every cell in the population not
    // belonging to the clone filter_clone, if negate_filter is true), this
    // will call the procedure proc.
    virtual void map(function<void (Cell &)> proc, 
                      const Clone * filter_clone = NULL,
                      bool negate_filter = false) = 0;
    // This is like map except the order of the cells is guaranteed to be
    // random.
    virtual void map_in_random_order(function<void (Cell &)> proc,
                                      const Clone * filter_clone = NULL,
                                      bool negate_filter = false) = 0;
    // This is like map except const.
    virtual void const_map(function<void (const Cell &)> proc, 
                            const Clone * filter_clone = NULL,
                            bool negate_filter = false) const = 0;
    // This is like map_in_random_order except const.
    virtual void const_map_in_random_order(function<void (const Cell &)> proc,
                                    const Clone * filter_clone = NULL,
                                    bool negate_filter = false) const = 0;
    // For every cell in the population belonging to the clone filter_clone (or
    // every cell, if filter_clone is NULL; or every cell in the population not
    // belonging to the clone filter_clone, if negate_filter is true), this
    // will call the function func.  The first time, it will pass
    // initial_value_ptr as the first argument to func.  Subsequently it will
    // pass the return value from the previous call to func as the first
    // argument to the next call to func.  It returns the return value of the
    // last call to func.
    virtual void * fold(function<void * (void *, Cell &)> func, 
                        void * initial_value_ptr, 
                        const Clone * filter_clone = NULL,
                        bool negate_filter = false) = 0;
    // This is like fold except the order of the cells is guaranteed to be
    // random.
    virtual void * fold_in_random_order(function<void * (void *, Cell &)> func,
                                        void * initial_value_ptr,
                                        const Clone * filter_clone = NULL,
                                        bool negate_filter = false) = 0;
    // This is like fold except const, as far as the population of cells is
    // concerned.
    virtual void * const_fold(function<void * (void *, const Cell &)> func, 
                              void * initial_value_ptr, 
                              const Clone * filter_clone = NULL,
                              bool negate_filter = false) const = 0;
    // This is like fold_in_random_order except const, as far as the population
    // of cells is concerned.
    virtual void * const_fold_in_random_order(
                                  function<void * (void *, const Cell &)> func,
                                  void * initial_value_ptr,
                                  const Clone * filter_clone = NULL,
                                  bool negate_filter = false) const = 0;
    // For each of the neighbors of the argument cell belonging to the clone
    // filter_neighbor_clone (or each of the neighbors of the argument cell, if
    // filter_neighbor_clone is NULL; or each of the neighbors of the argument
    // cell not belonging to the clone filter_clone, if negate_filter is true),
    // this will call the procedure proc with cell as the first argument and a
    // reference to the neighbor as the second argument.
    virtual void map_neighbors(function<void (Cell &, Cell &)> proc,
                                Cell &cell,
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false) = 0;
    // This is like map_neighbors except the order of the neighbors is
    // guaranteed to be random.
    virtual void map_neighbors_in_random_order(
                                function<void (Cell &, Cell &)> proc,
                                Cell &cell,
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false) = 0;
    // This is like map_neighbors, except const.
    virtual void const_map_neighbors(
                          function<void (const Cell &, const Cell &)> proc,
                          const Cell &cell,
                          const Clone * filter_neighbor_clone = NULL,
                          bool negate_filter = false) const = 0;
    // This is like map_neighbors_in_random_order, except const.
    virtual void const_map_neighbors_in_random_order(
                          function<void (const Cell &, const Cell &)> proc,
                          const Cell & cell_ptr,
                          const Clone * filter_neighbor_clone = NULL,
                          bool negate_filter = false) const = 0;
    // For each of the neighbors of the argument cell belonging to the clone
    // filter_neighbor_clone (or each of the neighbors of the argument cell, if
    // filter_neighbor_clone is NULL; or each of the neighbors of the argument
    // cell not belonging to filter_clone, if negate_filter is true), this will
    // call the function func with cell as the second argument and a reference
    // to the neighbor as the third argument.  The first time, it will pass
    // initial_value_ptr as the first argument to func.  Subsequently it will
    // pass the return value from the previous call to func as the first
    // argument to the next call to func.  It returns the return value of the
    // last call to func.
    virtual void * fold_neighbors(
                                function<void * (void *, Cell &, Cell &)> func, 
                                Cell &cell,
                                void * initial_value_ptr, 
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false) = 0;
    // This is like fold_neighbors except the order of the neighbors is
    // guaranteed to be random.
    virtual void * fold_neighbors_in_random_order(
                                function<void * (void *, Cell &, Cell &)> func, 
                                Cell &cell,
                                void * initial_value_ptr, 
                                const Clone * filter_neighbor_clone = NULL,
                                bool negate_filter = false) = 0;
    // This is like fold_neighbors except const, as far as the population of
    // cells is concerned.
    virtual void * const_fold_neighbors(
                    function<void * (void *, const Cell &, const Cell &)> func, 
                    const Cell &cell,
                    void * initial_value_ptr, 
                    const Clone * filter_neighbor_clone = NULL,
                    bool negate_filter = false) const = 0;
    // This is like fold_neighbors_in_random_order except const, as far as the
    // population of cells is concerned.
    virtual void * const_fold_neighbors_in_random_order(
                    function<void * (void *, const Cell &, const Cell &)> func, 
                    const Cell &cell,
                    void * initial_value_ptr, 
                    const Clone * filter_neighbor_clone = NULL,
                    bool negate_filter = false) const = 0;

    /*
     * End of methods to act globally or locally on the population of cells
     */

    virtual bool are_neighbors(const Cell &cell0, const Cell &cell1) const = 0;
    // This returns 0.0 if the cells are not neighbors.
    virtual float get_size_of_interface_between_cells(
                                                  const Cell &cell0, 
                                                  const Cell &cell1) const = 0;
    virtual bool is_adjacent_to_clone(const Cell &cell, 
                                      const Clone &clone) const;
    virtual int get_num_cells(Clone *filter_clone = NULL) const;
    virtual float get_volume(Clone *filter_clone = NULL) const;
    virtual int get_num_cells_of_clone_adjacent_to_clone(const Clone &clone,
                                      const Clone &neighbor_clone) const;
    virtual float get_size_of_boundary_of_clone(const Clone &clone) const;
    virtual float get_median_distance_from_clone_to_clone(const Clone &clone,
                                    const Clone &neighbor_clone) const = 0;
    virtual float get_max_fitness_ever(void) const { return max_fitness_ever; };
    virtual void kill(Cell &cell) {
      cell.kill();
    };
    virtual void update_all_fitnesses(void);
    // This returns NULL if there is no space to replicate, or a description of
    // the available space if there is.
    virtual void *check_for_space_to_replicate(Cell &cell) = 0;
    // This replicates cell into the specified space and records the
    // replication in replication_record.
    virtual void replicate(Cell &cell, void *space_specification,
                          ReplicationRecord &replication_record) = 0;
    function<void(Cell &, Cell &)> &get_neighbor_affect_cell_func(void) {
      return neighbor_affect_cell_func;
    }
    const std::set<float> &get_survival_probabilities(void) const { 
      return survival_probabilities; 
    };
    int get_num_cells_of_survival_probability(float survival_probability) 
        const {
      std::map<float, std::vector<Cell *> *>::const_iterator iter;
      iter = cells_of_survival_probability.find(survival_probability);
      if (iter != cells_of_survival_probability.end()) {
        return iter->second->size();
      }
      return 0;
    }
    Cell &get_cell_of_survival_probability(float survival_probability, int i) {
      return *(cells_of_survival_probability[survival_probability]->at(i));
    }
  protected:
    virtual void set_max_fitness_ever(float fitness) {
      max_fitness_ever = fitness;
    };
    void update_fitness(Cell &cell);
    friend void *update_fitness_func(void *data, Cell &cell);
  private:
    int initial_width, initial_height, initial_depth;
    float max_fitness_ever;
    // The second argument is the neighbor which is affecting the first
    // argument.
    function<void(Cell &, Cell &)> neighbor_affect_cell_func;
    void clear_cells_of_survival_probability(void);
    std::set<float> survival_probabilities;
    std::map<float, std::vector<Cell *> *> cells_of_survival_probability;
    // Only update_all_fitnesses() should call this.
};

#endif
