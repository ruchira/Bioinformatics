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
#include "distinct_hues.h"
#include "hsv_to_rgb.h"
#include <math.h>
#include <cassert>

VisualizeHexPopulation::VisualizeHexPopulation(int width, int height,
                                          std::vector<Clone *> &clone_ptrs,
                                  HexagonRendering &a_hexagon_rendering,
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

  fitness_increment = get_max_fitness_ever() / 255;

  strip_width_in_pixels = get_left_coord_in_pixels_of_cell_at(
                                    get_initial_width(), 0);
  frame_width_in_pixels = get_left_coord_in_pixels_of_cell_at(
                                   get_initial_width(), get_initial_height());
  frame_height_in_pixels = get_top_coord_in_pixels_of_cell_at(
                                   get_initial_width(), get_initial_height()-1)
                          + 2 * hexagon_rendering.get_side();

  int status = set_gfx_mode(GFX_AUTODETECT_WINDOWED, 640, 480,
                            frame_width_in_pixels, frame_height_in_pixels);
  status = install_keyboard();
  // Fill with black
  clear_bitmap(screen);
  for (int clone_index = 0; clone_index < clone_ptrs.size(); ++clone_index) {
    index_of_clone[clone_ptrs[clone_index]] = clone_index;
  }
  initialize_colors(clone_ptrs.size());
  frame = create_bitmap(frame_width_in_pixels, frame_height_in_pixels);
  clear_bitmap(frame);
  hexagon_rendering.initialize();
}

VisualizeHexPopulation::~VisualizeHexPopulation() {
  destroy_bitmap(frame);
}

void VisualizeHexPopulation::initialize_colors(int num_clones) {
  assert(num_clones <= 10);

  // The last entry in the palette is black, to represent a dead cell.  
  palette[HexagonRendering::black].r = 0;
  palette[HexagonRendering::black].g = 0;
  palette[HexagonRendering::black].b = 0;

  const int *distinct_hue_degrees = get_distinct_hue_degrees(num_clones);

  // Palette index 0 represents the mask color, so it is not actually used.
  palette[0].r = 0;
  palette[0].g = 0;
  palette[0].b = 0;
  int current_palette_index;
  // Leaving out the mask color 0 and the last color black, there are PAL_SIZE
  // colors in the palette; these are partitioned among the hues.  We've
  // already put the maximum level, with value = 1.0, in the palette.
  int num_levels = std::min(63, (PAL_SIZE - 2) / num_clones - 1);
  float value_increment = 1.0 / num_levels;
  int level;
  RGB rgb, rgb_increment;
  current_palette_index = num_clones + 1;
  for (int hue_num = 0; hue_num < num_clones; ++hue_num) {
    rgb.r = 0;
    rgb.g = 0;
    rgb.b = 0;
    hsv_to_rgb(distinct_hue_degrees[hue_num], 1.0, value_increment,
                rgb_increment);
    for (level = 1; level < num_levels; ++level) {
      rgb.r += rgb_increment.r;
      rgb.g += rgb_increment.g;
      rgb.b += rgb_increment.b;
      palette[current_palette_index].r = rgb.r;
      palette[current_palette_index].g = rgb.g;
      palette[current_palette_index].b = rgb.b;
      ++current_palette_index;
    }
    rgb.r += rgb_increment.r;
    rgb.g += rgb_increment.g;
    rgb.b += rgb_increment.b;
    palette[hue_num+1].r = std::min(63, (int)rgb.r);
    palette[hue_num+1].g = std::min(63, (int)rgb.g);
    palette[hue_num+1].b = std::min(63, (int)rgb.b);
  }
  for (; current_palette_index < PAL_SIZE - 1; ++current_palette_index) {
    palette[current_palette_index].r = 0;
    palette[current_palette_index].g = 0;
    palette[current_palette_index].b = 0;
  }
  set_palette(palette);
  create_light_table(&light_table, palette, 0, 0, 0, NULL);
  color_map = &light_table;
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
  int light, clone_number;
  if (cell.is_alive()) {
    light = ceil(cell.get_fitness() / fitness_increment);
    clone_number = 1 + index_of_clone[cell.get_clone_ptr()];
  } else {
    light = 0;
    clone_number = 0;
  }
  hexagon_rendering.render(frame, x, y, light, clone_number);
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
    hexagon_rendering.render(frame, x + strip_width_in_pixels, y, 255, 
                              clone_number);
  }
  if (x + hexagon_rendering.get_width() - strip_width_in_pixels > 0) {
    // The hexagon is in region C; render it again in region A
    hexagon_rendering.render(frame, x - strip_width_in_pixels, y, light,
                              clone_number);
  }
}

void VisualizeHexPopulation::envivify(void) {
  HexPopulation::envivify();
  fold(color_hex_cell_func, this);
}
