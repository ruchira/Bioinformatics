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
#include "hsv_to_rgb.h"
#include "distinct_hues.h"
#include <stdlib.h>
#include <cassert>
#include <iostream>
#include <sstream>

HexagonRendering::HexagonRendering(int new_side, int new_width, 
                                    int new_num_hues) : 
  side(new_side), width(new_width), num_hues(new_num_hues) {
  assert(new_side % 2 == 0);
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

}

void HexagonRendering::initialize(void) {
  specify_diagonal_line();

  // Initialize the palette to black. (Any unused colors will remain black.)
  set_palette(black_palette);


  // The last entry in the palette is black, to represent a dead cell.  Since
  // the whole palette is black already, we can skip initializing this one.

  assert(num_hues <= 10);

  RGB rgb;
  const int *distinct_hue_degrees = get_distinct_hue_degrees(num_hues);

  sprites = (BITMAP **)malloc((num_hues + 1) * sizeof(BITMAP *));
  int sprite_index;
  // Palette index 0 represents the mask color, so we skip it.
  int current_palette_index = 1;
  for (sprite_index = 0; sprite_index <= num_hues; ++sprite_index) {
    sprites[sprite_index] = create_bitmap(width, 2 * side);
    // This clears the bitmap to the mask color
    clear_bitmap(sprites[sprite_index]);
    if (sprite_index == 0) {
      // Fill the hexagon with black
      fill_hexagon(sprite_index, black);
    } else {
      hsv_to_rgb(distinct_hue_degrees[sprite_index - 1], 1.0, 1.0, rgb);
      set_color(current_palette_index, &rgb);
      fill_hexagon(sprite_index, current_palette_index);
      ++current_palette_index;
    }
  }
  PALETTE palette;
  get_palette(palette);
  for (sprite_index = 0; sprite_index <= num_hues; ++sprite_index) {
    std::ostringstream os;
    os << "sprite_" << sprite_index << ".bmp";
    save_bmp(os.str().c_str(), sprites[sprite_index], palette);
  }
  create_light_table(&light_table, palette, 0, 0, 0, NULL);
}

HexagonRendering::~HexagonRendering() {
  for (int sprite_index = 0; sprite_index <= num_hues; ++sprite_index) {
    destroy_bitmap(sprites[sprite_index]);
  }
  free(diagonal_line_specification);
}

void HexagonRendering::fill_hexagon(int sprite_index, int color) {
  int row, upper_row, lower_row, min_col, max_col;
  // The vertical sides are the leftmost and rightmost columns, going from
  // vertical side/2 through side/2 + side - 1.  The leftmost and rightmost
  // columns have coordinates 0 and width - 1, so the *interior* of the hexagon
  // has horizontal coordinates from 1 to width - 2 inclusive.
  min_col = 1;
  max_col = width - 2;
  for (row = side/2; row < side/2 + side; ++row) {
    hline(sprites[sprite_index], min_col, row, max_col, color);
  }
  int i;
  for (i = 0, upper_row = side/2 - 1, lower_row = side/2 + side;
      i < side/2 - 1; 
      i++, upper_row--, lower_row++) {
    min_col += diagonal_line_specification[i];
    max_col -= diagonal_line_specification[i];
    hline(sprites[sprite_index], min_col, upper_row, max_col, color);
    hline(sprites[sprite_index], min_col, lower_row, max_col, color);
  }
}
