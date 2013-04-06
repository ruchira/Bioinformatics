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
  public:
    HexCell(): horiz_coord(INT_MIN), diag_coord(INT_MIN) {};
    virtual ~HexCell() {};
    int get_horiz_coord(void) const { return horiz_coord; };
    void set_horiz_coord(int horiz) { horiz_coord = horiz; };
    int get_diag_coord(void) const { return diag_coord; };
    void set_diag_coord(int diag) { diag_coord = diag; };
    virtual float get_volume(void) const { return volume; };
  private:
    static const float volume;
    int horiz_coord, diag_coord;
};

#endif
