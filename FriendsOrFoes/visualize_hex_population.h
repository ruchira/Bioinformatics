// visualize_hex_population.h
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
#ifndef VISUALIZE_HEX_POPULATION_H
#define VISUALIZE_HEX_POPULATION_H
#include "hex_population.h"
#include "hexagon_rendering.h"
#include <allegro.h>
#include <vector>
#include <map>
#include <iostream>

class VisualizeHexPopulation : public HexPopulation {
	public:
		VisualizeHexPopulation(int width, int height, 
                            std::vector<Clone *> &clone_ptrs,
                            HexagonRendering &a_hexagon_rendering,
                            float max_fitness_ever_seen);
		virtual ~VisualizeHexPopulation();
    // This will write the current frame to the named output file.
		void visualize(const char *output_filename) {
      acquire_bitmap(frame);
      int status = save_bitmap(output_filename, frame, palette);
      release_bitmap(frame);
      draw_sprite(screen, frame, 0, 0);
    };
    virtual void envivify(void);
    virtual void kill(Cell &cell) { 
      HexPopulation::kill(cell);
      color_cell(cell);
    };
    virtual void update_fitness(Cell &cell) {
      HexPopulation::update_fitness(cell);
      color_cell(cell);
    };
    virtual void replicate(Cell &cell, void *space_specification,
                          ReplicationRecord &replication_record) {
      HexPopulation::replicate(cell, space_specification, replication_record);
      color_cell(*replication_record.get_daughter0_ptr());
      color_cell(*replication_record.get_daughter1_ptr());
    };
  protected:
    virtual void initialize_colors(int num_clones);
    int get_left_coord_in_pixels_of_cell_at(int horiz_coord, int diag_coord) {
      // The rightmost column of one hexagon is the leftmost column of the
      // next, so each hexagon effectively takes up only
      // hexagon_rendering.get_width() - 1 pixels in width.
      return (horiz_coord * (hexagon_rendering.get_width() - 1) +
            diag_coord 
            * hexagon_rendering.get_horizontal_component_of_diagonal_offset());
    }
    int get_top_coord_in_pixels_of_cell_at(int horiz_coord, int diag_coord) {
      return diag_coord 
            * hexagon_rendering.get_vertical_component_of_diagonal_offset();
    }
	private:
    PALETTE palette;
    int num_values;
    std::map<Clone *, int> index_of_clone;
    float fitness_increment;
    HexagonRendering &hexagon_rendering;
    // Each clone has a distinct hue; cells are colored with decreasing
    // value (as in hue, saturation, value)  as they become less fit.
    // These take const arguments because the cell instance is the *model* of
    // the cell, whereas these actually just color the *view* of the cell in
    // the rendered frame.
    void color_hex_cell(const HexCell &cell);
    void color_cell(const Cell &cell) {
      color_hex_cell( *((const HexCell *)(&cell)) );
    }
    friend void * color_hex_cell_func(void *data, Cell &cell);
		BITMAP *frame;
    int frame_width_in_pixels;
    int frame_height_in_pixels;
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
    int strip_width_in_pixels;
    COLOR_MAP light_table;
};

#endif
