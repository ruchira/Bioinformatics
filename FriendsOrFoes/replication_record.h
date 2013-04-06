// replication_record.h
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
#ifndef REPLICATION_RECORD_H
#define REPLICATION_RECORD_H
#include "cell.h"

class ReplicationRecord {
  public:
    ReplicationRecord() : mother_ptr(NULL), daughter0_ptr(NULL),
                                            daughter1_ptr(NULL) {};
    virtual ~ReplicationRecord() {
      if (mother_ptr != NULL) {
        delete mother_ptr;
      }
    }
    // This assigns a new object of a concrete subclass of Cell to mother_ptr.
    virtual void initialize(void) = 0;
    const Cell *get_mother_ptr(void) const { return mother_ptr; };
    const Cell *get_daughter0_ptr(void) const { return daughter0_ptr; };
    const Cell *get_daughter1_ptr(void) const { return daughter1_ptr; };
    // This copies the data of cell_ptr as a concrete subclass of Cell, to
    // the object pointed to by mother_ptr.
    virtual void record_as_mother(const Cell *cell_ptr) = 0;
    void set_daughter0_ptr(Cell *cell_ptr) { daughter0_ptr = cell_ptr; };
    void set_daughter1_ptr(Cell *cell_ptr) { daughter1_ptr = cell_ptr; };
  protected:
    void set_mother_ptr(Cell *cell_ptr) { mother_ptr = cell_ptr; };
  private:
    // The replication record owns the Cell instance pointed to by mother_ptr
    // and is responsible for deleting it.
    Cell *mother_ptr;
    // The replication record does not own the Cell instances pointed to by
    // daughter0_ptr and daughter1_ptr, and should not change them.
    Cell *daughter0_ptr;
    Cell *daughter1_ptr;
};

#endif
