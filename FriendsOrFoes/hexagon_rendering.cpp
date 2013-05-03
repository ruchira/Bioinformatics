// hexagon_rendering.cpp
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
#include "hexagon_rendering.h"
#include <stdlib.h>
#include <cassert>

HexagonRendering::HexagonRendering(int new_side, int new_width) {
  if (new_side % 2 == 0) {
    side = new_side;
    assert(new_width % 2 == 1 && ((new_width - 1) / 2) % 2 == (side / 2) % 2);
    width = new_width;
    /*
    |   |   |
    | a |   |
     \ / \ /
      |   |
      |   |
      | b |
     / \ / \
    |   |   |
    |   |   |
    The diagonal offset is between the point a, the upper left corner of the
    rectangular region containing the middle hexagon, and the point b, the
    corresponding point in the rectangular region containing the lower right
    hexagon.
    */
    diagonal_offset.first = (width - 1) / 2;
    diagonal_offset.second = side/2 + side - 1;
    diagonal_line_specification = (int *)malloc((side / 2 - 1) * sizeof(int));
    specify_diagonal_line();
    sprite = create_bitmap(width - 1, 2 * side - 1);
    // This clears the bitmap with color 0, the mask color (which is white)
    clear_bitmap(sprite);
    fill_hexagon(black);
  }
}

HexagonRendering::~HexagonRendering() {
  destroy_bitmap(sprite);
  free(diagonal_line_specification);
}

void HexagonRendering::fill_hexagon(int color) {
  int row, upper_row, lower_row, min_col, max_col;
  // The vertical sides are the leftmost and rightmost columns, going from
  // vertical side/2 through side/2 + side - 1.  The leftmost and rightmost
  // columns have coordinates 0 and width - 1, so the *interior* of the hexagon
  // has horizontal coordinates from 1 to width - 2 inclusive.
  min_col = 1;
  max_col = width - 2;
  for (row = side/2; row < side/2 + side; ++row) {
    hline(sprite, min_col, row, max_col, color);
  }
  int i;
  for (i = 0, upper_row = side/2 - 1, lower_row = side/2 + side;
      i < side/2 - 1; 
      i++, upper_row--, lower_row++) {
    min_col += diagonal_line_specification[i];
    max_col -= diagonal_line_specification[i];
    hline(sprite, min_col, upper_row, max_col, color);
    hline(sprite, min_col, lower_row, max_col, color);
  }
}
