// cell_cycle.h
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
#ifndef CELL_CYCLE_H
#define CELL_CYCLE_H
#include "probability.h"
#include "population.h"
#include <vector>

using std::tr1::function;

class CellCycle{
  public:
    CellCycle(Population &a_population, 
          function<ReplicationRecord *(void)> new_replication_record_maker): 
      population(a_population), total_num_killed_cells(0),
      total_num_cells_that_replicated(0),
      make_new_replication_record(new_replication_record_maker) {};
    virtual ~CellCycle();
    virtual void run(void);
    const Population &get_population(void) const { return population; };
    const vector<const Cell *> *get_killed_cells_of_clone(const Clone &clone);
    // For every ReplicationRecord from the last run, this will call the
    // procedure proc.
    virtual void const_map_replication_records(
                        function<void (const ReplicationRecord &)> proc) const;
    // For every ReplicationRecord from the last run, this will call the
    // function func.  The first time, it will pass initial_value_ptr as the
    // second argument to func.  Subsequently it will pass the return value
    // from the previous call to func as the second argument to the next call
    // to func.  It returns the return value of the last call to func.
    virtual void * const_fold_replication_records(
                      function<void *(void *, const ReplicationRecord &)> func, 
                      void *initial_value_ptr) const;
    int get_total_num_killed_cells(void) const { 
      return total_num_killed_cells;
    };
    int get_total_num_cells_that_replicated(void) const {
      return total_num_cells_that_replicated;
    };
    virtual void check_survival(Cell &cell);
    virtual void check_reproduction(Cell &cell);
    virtual void global_check_reproduction(void);
  private:
    Population &population;
    function<ReplicationRecord *(void)> make_new_replication_record;
    map<const Clone *, vector<const Cell *> * > killed_cells_of_clone;
    map<const Clone *, vector<ReplicationRecord *> * > 
        replication_records_of_clone;
    // Each ReplicationRecord instance controls a block of memory
    // corresponding to a concrete subclass of Cell.  To avoid continually
    // allocating and deallocating this memory at every run, we don't actually
    // destroy the ReplicationRecords at the end of a run, and we keep these
    // instances in the vector.  This means the size of the vector does not
    // accurately reflect the number of replication records from the last run.
    // Instead, we keep track of that number here.
    map<const Clone *, int > num_replication_records_of_clone;
    int total_num_killed_cells;
    int total_num_cells_that_replicated;
    void clear_killed_cells_of_clone();
    void clear_num_replication_records_of_clone();
    ReplicationRecord &get_next_replication_record_ptr_of_clone(
                                                      const Clone &clone);
};

#endif
