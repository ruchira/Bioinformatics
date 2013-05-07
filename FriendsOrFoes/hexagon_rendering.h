// hexagon_rendering.h
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
#ifndef HEXAGON_RENDERING_H
#define HEXAGON_RENDERING_H

#include <utility>
#include <allegro.h>
#include <iostream>

class HexagonRendering {
  public:
    void render(BITMAP *destination, int x, int y, int light, 
                int sprite_number) const {
      draw_lit_sprite(destination, sprites[sprite_number], x, y, light);
    }
    int get_side(void) const { return side; };
    int get_width(void) const { return width; };
    int get_height(void) const { return 2 * side; };
    // Consider the x-axis to go from left to right and the y-axis to go from
    // top to bottom.  These methods return the coordinates of the top
    // left-corner of the next nearest hexagon, downward and to the right of
    // this one.
    int get_horizontal_component_of_diagonal_offset(void) const { 
      return diagonal_offset.first; 
    };
    int get_vertical_component_of_diagonal_offset(void) const {
      return diagonal_offset.second;
    };
    // The last color in the palette is black
    static const int black = PAL_SIZE - 1;
  protected:
    // The hexagon will occupy a rectangular region, 
    // width pixels wide by 2 * side pixels high.  
    // width / side should be close to sqrt(3).
    HexagonRendering(int new_side, int new_width, int num_hues);
    virtual ~HexagonRendering();
    virtual void initialize(void);
    // The leftmost and rightmost columns contain the vertical
    // sides of the hexagon.  Thus, here pixels with vertical coordinates from
    // side / 2 to side / 2 + side - 1 constitute the boundary.
    // The pixel with coordinates ((width - 1) / 2, 0) is the topmost
    // boundary pixel of the hexagon (and also, the bottom pixel of the
    // vertical edge of the nearest neighboring hexagons, to the upper left and
    // upper right).  The top left diagonal edge of the hexagon, then, goes
    // from (0, side / 2) up to ( (width - 1) / 2, 0).  It has to start
    // diagonally from either end (otherwise the vertical side would be longer
    // than it is), so the edge goes from (1, side / 2 - 1) up to 
    // ( (width - 1) / 2 - 1, 1 ).  It has to traverse a horizontal distance of
    // (width - 1) / 2 - 1 pixels and a vertical distance of side / 2 - 1
    // pixels.  It will do this by going through a series of horizontal runs of
    // varying lengths, adding
    // up to total horizontal distance (width - 1) / 2 - 1, after each of which
    // it steps up by 1 pixel.  The diagonal_line_specification specifies the
    // lengths of these (side / 2 - 1) horizontal runs.
    int *diagonal_line_specification; 
    // The concrete subclass will specify the rendering of the hexagon by
    // filling in this method, which will initialize the
    // diagonal_line_specification array.  The memory for the array is managed
    // by this HexagonRendering superclass, the subclass need only fill in the
    // values.  For symmetry the values should be palindromic.
    virtual void specify_diagonal_line(void) = 0;
  private:
    int side; // The number of pixels in the unit length.  
    int width;  // The number of pixels wide.
    int num_hues; // The number of different hues (besides black) to use

    // The offset of the upper left corner of the rectangular region
    // corresponding to the next lower right hexagonal neighbor of this
    // hexagon.
    std::pair<int, int> diagonal_offset;

    // These are the actual sprites, with the interior pixels of the rendered
    // hexagon colored and all other pixels at color 0 (indicating
    // transparency).
    // Note that the rightmost column and bottommost row are omitted, since
    // only black boundary pixels are included there from this hexagon.
    // So each sprite is actually width - 1 pixels wide and 2 * side - 1 pixels
    // high.
    // There are num_hues + 1 of these (the 0th one is black)
    BITMAP **sprites;
    COLOR_MAP light_table;
    // This fills the interior of the hexagon with color.
    void fill_hexagon(int sprite_index, int color);
};

#endif
