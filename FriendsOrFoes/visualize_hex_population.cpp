// visualize_hex_population.cpp
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
#include "visualize_hex_population.h"
#include <math.h>
#include <cassert>
#include "hsv_to_rgb.h"
#include "distinct_hues.h"

void VisualizeHexPopulation::initialize_colors(
                                            std::vector<Clone *> &clone_ptrs) {
  int num_clones = clone_ptrs.size();
  assert(num_clones <= 10);
  // Initialize the palette to black. (Any unused colors will remain black.)
  set_palette(black_palette);

  RGB rgb, rgb_increment;
  // The zeroth entry in the palette is the mask (transparent) color.
  // We'll make this white.
  rgb.r = 0xFF;
  rgb.g = 0xFF;
  rgb.b = 0xFF;
  set_color(0, &rgb);
  // The last entry in the palette is black, to represent a dead cell.  Since
  // the whole palette is black already, we can skip initializing this one.

  // Divide up the PAL_SIZE - 2 colors remaining evenly among clones.
  // Each clone will have a distinct hue.  The fitness level of a cell of a
  // particular clone will be depicted by its value (as in, hue, saturation,
  // value).
  num_values = (PAL_SIZE - 2) / num_clones;
  float value_increment = 1.0 / num_values;
  fitness_increment = get_max_fitness_ever() * value_increment;

  int current_palette_index = 1;
  const int *distinct_hue_degrees = get_distinct_hue_degrees(num_clones);
  int hue;

  int level;
  float current_value;
  for (int clone_num = 0; clone_num < num_clones; ++clone_num) {
    index_of_clone[clone_ptrs[clone_num]] = clone_num;
    hue = distinct_hue_degrees[clone_num];
    rgb.r = 0;
    rgb.g = 0;
    rgb.b = 0;
    hsv_to_rgb(hue, 1.0, value_increment, rgb_increment);
    for (level = 1; level <= num_values; ++level) {
      if (rgb.r + rgb_increment.r > 255) {
        rgb.r = 255;
      } else {
        rgb.r += rgb_increment.r;
      }
      if (rgb.g + rgb_increment.g > 255) {
        rgb.g = 255;
      } else {
        rgb.g += rgb_increment.g;
      }
      if (rgb.b + rgb_increment.b > 255) {
        rgb.b = 255;
      } else {
        rgb.b += rgb_increment.b;
      }
      set_color(current_palette_index, &rgb);
      ++current_palette_index;
    }
  }
}

VisualizeHexPopulation::VisualizeHexPopulation(int width, int height,
                                          std::vector<Clone *> &clone_ptrs,
                                  const HexagonRendering &a_hexagon_rendering,
                                  float max_fitness_ever_seen) :
	HexPopulation(width, height), hexagon_rendering(a_hexagon_rendering) {
  set_max_fitness_ever(max_fitness_ever_seen);
	// The hexagons have sides of unit length.  This unit length is to be
	// represented by num_pixels_in_unit_length pixels.  Each hexagon has
	// two vertical sides.  
	// In HexPopulation we have the following diagram:
	/*
     5---4        ^ increasing row
     /\ /\       /
    6--o--3     /
     \/ \/     /
     1---2     ----> increasing column
    The horizontal dimension wraps around, the diagonal dimension doesn't, thus 
    making the grid topologically a tube like the esophagus
  */
	// Here o, 1, 2, ..., 6 are centers of neighboring hexagons.  Thus, the sides
	// of the hexagon centered on o are perpendicular bisectors of the lines of
	// the above diagram radiating from o:
	/*
		  A
		 / \
    |   |
		| o |
    |   |
		 \ /
		  B
	*/
	// The horizontal distance from the center of this hexagon to one of the
	// vertical sides is the altitude of an equilateral triangle of side one,
	// i.e., it is sqrt(3)/2.  So the width of a single hexagon is sqrt(3).
	// The vertical distance from the center of this hexagon to the vertex A is
	// the length of an equilateral triangle of side one, and so is the distance
	// from the center to vertex B.  So the height of a single hexagon is 2.

  initialize_colors(clone_ptrs);

  strip_width_in_pixels = get_left_coord_in_pixels_of_cell_at(
                                    get_initial_width(), 0);
  frame_width_in_pixels = get_left_coord_in_pixels_of_cell_at(
                                   get_initial_width(), get_initial_height());
  int height_in_pixels = get_top_coord_in_pixels_of_cell_at(get_initial_width(),
                                                          get_initial_height());

  frame = create_bitmap(frame_width_in_pixels, height_in_pixels);
  // Fill with black
  rectfill(frame, 0, 0, frame_width_in_pixels - 1, height_in_pixels - 1,
            HexagonRendering::black);
}

VisualizeHexPopulation::~VisualizeHexPopulation() {
  destroy_bitmap(frame);
}

void * color_hex_cell_func(void *data, Cell &cell) {
  VisualizeHexPopulation *population_ptr = (VisualizeHexPopulation *)data;
  population_ptr->color_cell(cell);
  return data;
}

void VisualizeHexPopulation::color_hex_cell(const HexCell &cell) {
  int x = get_left_coord_in_pixels_of_cell_at(cell.get_horiz_coord(),
                                              cell.get_diag_coord());
  int y = get_top_coord_in_pixels_of_cell_at(cell.get_horiz_coord(),
                                              cell.get_diag_coord());
  int color;
  if (cell.is_alive()) {
    int i = index_of_clone[cell.get_clone_ptr()];
    int level = ceil(cell.get_fitness() / fitness_increment);
    color = level + i * num_values;
  } else {
    color = HexagonRendering::black;
  }
  hexagon_rendering.render(frame, x, y, color);
  // The hexagon coordinates run horizontally and diagonally within the
  // rectangular frame.  So the strip of hexagons is positioned like this
  // within the rectangular frame:
  /*
     -------------------
     |\ B:         :\ D|
     | \ :         : \ |
     |A \:         :C \|
     -------------------
   */
  // Since the strip wraps around horizontally, some of the hexagons need to
  // appear twice in the frame: the hexagons in region C of the strip appear
  // again in the frame in region A, and the hexagons in region B of the
  // strip appear again in the frame in region D.
  if (x + strip_width_in_pixels < frame_width_in_pixels) {
    // The hexagon is in region B; render it again in region D
    hexagon_rendering.render(frame, x + strip_width_in_pixels, y, color);
  }
  if (x + hexagon_rendering.get_width() - strip_width_in_pixels > 0) {
    // The hexagon is in region C; render it again in region A
    hexagon_rendering.render(frame, x - strip_width_in_pixels, y, color);
  }
}

void VisualizeHexPopulation::envivify(void) {
  HexPopulation::envivify();
  fold(color_hex_cell_func, this);
}
