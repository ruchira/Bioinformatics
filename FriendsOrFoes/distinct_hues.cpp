// distinct_hues.cpp
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
#include "distinct_hues.h"
#include <stdlib.h>

// These perceptually distinct hues were obtained from
// http://phrogz.net/css/distinct-colors.html
// fixing saturation and value at 1, and allowing hue to vary from 0 to 360
// degrees with fixed step size in degrees.  The website perturbs the even
// spacing of colors around the color wheel if necessary, using the
// CMC(I:c) algorithm to ensure visually distinct colors.  See
// http://en.wikipedia.org/wiki/Color_difference .
// Number of colors   Step Size in degrees
// 2                  179
// 3                  119
// 4                  91
// 5                  72
// 6                  59
// 7                  40
// 8                  37                
// 9                  28
// 10                 25
const int distinct_hue_degrees_1[1] = { 0 };
const int distinct_hue_degrees_2[2] = { 0, 179 };
const int distinct_hue_degrees_3[3] = { 0, 119, 238 };
const int distinct_hue_degrees_4[4] = { 0, 91, 182, 273 };
const int distinct_hue_degrees_5[5] = { 0, 72, 144, 216, 288 };
const int distinct_hue_degrees_6[6]
  = { 0, 59, 118, 177, 236, 295 };
const int distinct_hue_degrees_7[7]
  = { 0, 40, 120, 160, 200, 240, 320 };
const int distinct_hue_degrees_8[8] 
  = { 0, 37, 74, 148, 185, 222, 296, 333 };
const int distinct_hue_degrees_9[9] 
  = { 0, 28, 56, 84, 168, 196, 252, 308, 336 };
const int distinct_hue_degrees_10[10]
  = { 0, 25, 50, 75, 125, 175, 200, 275, 300, 325 };

const int *get_distinct_hue_degrees(int num_hues) {
  switch (num_hues) {
    case 1:
      return distinct_hue_degrees_1;
    case 2:
      return distinct_hue_degrees_2;
    case 3:
      return distinct_hue_degrees_3;
    case 4:
      return distinct_hue_degrees_4;
    case 5:
      return distinct_hue_degrees_5;
    case 6:
      return distinct_hue_degrees_6;
    case 7:
      return distinct_hue_degrees_7;
    case 8:
      return distinct_hue_degrees_8;
    case 9:
      return distinct_hue_degrees_9;
    case 10:
      return distinct_hue_degrees_10;
    default:
      return NULL;
  };
  return NULL;
}
