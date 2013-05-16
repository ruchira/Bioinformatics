// hex_cell.h
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
#ifndef HEX_CELL_H
#define HEX_CELL_H
#include "cell.h"
#include <limits.h>

class HexCell: public Cell {
  friend class HexPopulation;
  public:
    HexCell(): horiz_coord(INT_MIN), diag_coord(INT_MIN) {};
    virtual ~HexCell() {
      if (non_null_neighbor != NULL) {
        free(non_null_neighbor);
        non_null_neighbor = NULL;
      }
    };
    int get_horiz_coord(void) const { return horiz_coord; };
    void set_horiz_coord(int horiz) { horiz_coord = horiz; };
    int get_diag_coord(void) const { return diag_coord; };
    void set_diag_coord(int diag) { diag_coord = diag; };
    virtual float get_volume(void) const { return volume; };
    int get_num_neighbors(void) const { return num_neighbors; };
    // The neighbors are numbered in a fixed order:
    /* There are six possible neighbors of each cell:
     5---4        ^ increasing row
     /\ /\       /
    6--o--3     /
     \/ \/     /
     1---2     ----> increasing column
    The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
    making the grid topologically a tube like the esophagus
    */
    // These functions return NULL if the cell has no neighbor at the given
    // index (for example, if the cell is at the bottom of the tube, it has no
    // neighbor 1 or 2).
    HexCell *get_neighborptr(int i) { return neighbor[i-1]; };
    const HexCell *const_get_neighborptr(int i) const { return neighbor[i-1]; };
    // These functions return the ith non-null neighbor numbered from 1 to
    // num_neighbors.  For example, if the cell is at the bottom edge and thus
    // has 4 neighbors, namely neighbors 3, 4, 5, and 6 in the above diagram,
    // then: get_non_null_neighbor(1) will return neighbor 3, ...
    // get_non_null_neighbor(4) will return neighbor 6.  If the cell is at the
    // top edge and thus has 4 neighbors, namely neighbors 1, 2, 3, and 6 in
    // the above diagram, then get_non_null_neighbor(1) will return neighbor 1,
    // ... get_non_null_neighbor(3) will return neighbor 3, and
    // get_non_null_neighbor(4) will return neighbor 6.
    HexCell &get_non_null_neighbor(int i) {
      return (num_neighbors == 6 ? *neighbor[i-1] : *non_null_neighbor[i-1]);
    }
    const HexCell &const_get_non_null_neighbor(int i) const {
      return (num_neighbors == 6 ? *neighbor[i-1] : *non_null_neighbor[i-1]);
    }
  private:
    static const float volume;
    int horiz_coord, diag_coord;
    // The HexPopulation friend class is responsible for maintaining these.
    int num_neighbors;
    HexCell* neighbor[6];
    // If num_neighbors < 6, this will be non-null and will be an array of the
    // non-null-neighbors.  
    HexCell** non_null_neighbor;
};

#endif
