// cell_cycle.cpp
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
#include "cell_cycle.h"

void CellCycle::kill(Cell &cell) {
  std::map<const Clone *, std::vector<const Cell *> * >::const_iterator iter;
  iter = killed_cells_of_clone.find(cell.get_clone_ptr());
  if (iter == killed_cells_of_clone.end()) {
    killed_cells_of_clone[cell.get_clone_ptr()] 
      = new std::vector<const Cell *>;
  }
  killed_cells_of_clone[cell.get_clone_ptr()]->push_back(&cell);
  population.kill(cell);
  ++total_num_killed_cells;
}

void CellCycle::try_survival(Cell &cell) {
  if (cell.is_alive()) {
    double survival_probability;
    survival_probability = clamp_probability(
                            cell.get_clone_ptr()->get_survival_coefficient() 
                            * cell.get_fitness());
    if (survival_probability == 1.0) {
      return;
    }
    if (survival_probability == 0.0 
        || !Random::bernoulli_check(survival_probability)) {
      kill(cell);
    }
  }
}

void *try_survival_func(void *data, Cell &cell) {
  CellCycle *cell_cycle_ptr = (CellCycle *)data;
  cell_cycle_ptr->try_survival(cell);
  return data;
}

void CellCycle::global_try_survival(void) {
  const std::set<float>& survival_probabilities 
    = population.get_survival_probabilities();
  std::set<float>::const_iterator prob_iter;
  std::set<int>::const_iterator index_iter;
  int n, i, j;
  unsigned int num_surviving_cells, num_cells_to_kill;
  for (prob_iter = survival_probabilities.begin(); 
        prob_iter != survival_probabilities.end(); ++prob_iter) {
    if (*prob_iter < 1.0) {
      n = population.get_num_cells_of_survival_probability(*prob_iter);
      if (n > 0) {
        // Do a binomial check to see how many of these cells survive.
        num_surviving_cells = Random::binomial_check(*prob_iter, n);
        if (num_surviving_cells < n) {
          if (num_surviving_cells == 0) {
            for (i = 0; i < n; ++i) {
              kill(population.get_cell_of_survival_probability(*prob_iter, i));
            }
          } else {
            num_cells_to_kill = n - num_surviving_cells;
            indices_of_cells_to_kill.clear();
            // Use Floyd's algorithm to pick a uniform random subset of size
            // num_cells_to_kill of the cells with this survival probability.
            for (i = n - num_cells_to_kill; i < n; ++i) {
              j = Random::rand_int_between_inclusive(0, i);
              index_iter = indices_of_cells_to_kill.find(j);
              if (index_iter == indices_of_cells_to_kill.end()) {
                indices_of_cells_to_kill.insert(j);
              } else {
                indices_of_cells_to_kill.insert(i);
              }
            }
            for (index_iter = indices_of_cells_to_kill.begin();
                  index_iter != indices_of_cells_to_kill.end(); ++index_iter) {
              kill(population.get_cell_of_survival_probability(*prob_iter,
                                                              *index_iter));
            }
          }
        }
      }
    }
  }
}

void CellCycle::try_reproduction(Cell &cell) {
  void *space_specification = population.check_for_space_to_replicate(cell);
  if (space_specification != NULL) {
    double reproduction_probability;
    reproduction_probability 
      = clamp_probability(cell.get_clone_ptr()->get_reproduction_coefficient() 
                          * cell.get_fitness());
    if (Random::bernoulli_check(reproduction_probability)) {
      population.replicate(cell, space_specification,
                get_next_replication_record_ptr_of_clone(*cell.get_clone_ptr()));
      ++total_num_cells_that_replicated;
    }
  }
}

void *try_reproduction_func(void *data, Cell &cell) {
  CellCycle *cell_cycle_ptr = (CellCycle *)data;
  cell_cycle_ptr->try_reproduction(cell);
  return data;
}

void CellCycle::global_try_reproduction(void) {
  population.fold_in_random_order(try_reproduction_func, this);
}

void CellCycle::run(void) {
  clear_killed_cells_of_clone();
  clear_num_replication_records_of_clone();
  total_num_killed_cells = 0;
  total_num_cells_that_replicated = 0;
  // We update the fitness on all the cells before checking survival, so that
  // the aliveness check will refer consistently to the previous run.
  population.update_all_fitnesses();
  global_try_survival();
  global_try_reproduction();
}

const std::vector<const Cell *> *CellCycle::get_killed_cells_of_clone(
                                                          const Clone &clone) {
  std::map<const Clone *, std::vector<const Cell *> *>::const_iterator iter;
  iter = killed_cells_of_clone.find(&clone);
  if (iter == killed_cells_of_clone.end()) {
    return NULL;
  } else {
    return iter->second;
  }
}

void CellCycle::clear_killed_cells_of_clone(void) {
  std::map<const Clone *, std::vector<const Cell *> *>::const_iterator iter;
  for (iter = killed_cells_of_clone.begin(); 
      iter != killed_cells_of_clone.end(); ++iter) {
    iter->second->clear();
  }
}

void CellCycle::clear_num_replication_records_of_clone(void) {
  std::map<const Clone *, int>::const_iterator iter;
  for (iter = num_replication_records_of_clone.begin(); 
      iter != num_replication_records_of_clone.end(); ++iter) {
    num_replication_records_of_clone[iter->first] = 0;
  }
}

void CellCycle::const_map_replication_records(
                        function<void (const ReplicationRecord &)> proc) const {
  std::map<const Clone *, int>::const_iterator num_iter;
  int i;
  for (num_iter = num_replication_records_of_clone.begin();
        num_iter != num_replication_records_of_clone.end();
        ++num_iter) {
    for (i = 0; i < num_iter->second; ++i) {
      proc(*replication_records_of_clone.at(num_iter->first)->at(i));
    }
  }
}

void *CellCycle::const_fold_replication_records(
                      function<void *(void *, const ReplicationRecord &)> func, 
                      void *initial_value_ptr) const {
  std::map<const Clone *, int>::const_iterator num_iter;
  int i;
  void *value_ptr = initial_value_ptr;
  for (num_iter = num_replication_records_of_clone.begin();
        num_iter != num_replication_records_of_clone.end();
        ++num_iter) {
    for (i = 0; i < num_iter->second; ++i) {
      value_ptr = func(value_ptr,
                       *replication_records_of_clone.at(num_iter->first)->at(i));
    }
  }
}

ReplicationRecord &CellCycle::get_next_replication_record_ptr_of_clone(
                                                        const Clone &clone) {
  std::map<const Clone *, std::vector<ReplicationRecord *> *>::const_iterator iter;
  iter = replication_records_of_clone.find(&clone);
  if (iter == replication_records_of_clone.end()) {
    replication_records_of_clone[&clone] 
        = new std::vector<ReplicationRecord *>;
  }
  std::map<const Clone *, int>::const_iterator num_iter;
  num_iter = num_replication_records_of_clone.find(&clone);
  if (num_iter == num_replication_records_of_clone.end()) {
    num_replication_records_of_clone[&clone] = 0;
  }
  ReplicationRecord *replication_record_ptr;
  if (num_replication_records_of_clone[&clone] <
      replication_records_of_clone[&clone]->size()) {
    replication_record_ptr 
      = replication_records_of_clone[&clone]->at(
                        num_replication_records_of_clone[&clone]);
  } else {
    replication_record_ptr = make_new_replication_record();
    replication_records_of_clone[&clone]->push_back(
                                                  replication_record_ptr);
  }
  num_replication_records_of_clone[&clone]++;
  return *replication_record_ptr;
}

CellCycle::~CellCycle() {
  std::map<const Clone *, std::vector<const Cell *> *>::const_iterator iter1;
  for (iter1 = killed_cells_of_clone.begin();
      iter1 != killed_cells_of_clone.end(); ++iter1) {
    delete iter1->second;
  }
  std::map<const Clone *, std::vector<ReplicationRecord *> *>::const_iterator iter2;
  for (iter2 = replication_records_of_clone.begin();
      iter2 != replication_records_of_clone.end(); ++iter2) {
    for (int i = 0; i < iter2->second->size(); ++i) {
      delete iter2->second->at(i);
    }
    delete iter2->second;
  }
}
