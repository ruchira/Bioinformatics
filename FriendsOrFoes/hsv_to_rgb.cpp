// hsv_to_rgb.cpp
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

#include "hsv_to_rgb.h"
#include <math.h>

void hsv_to_rgb(int hue_in_degrees, float saturation, float value, 
                        RGB &rgb) {
  // This implements the algorithm given on the Wikipedia page
  // http://en.wikipedia.org/wiki/HSL_and_HSV

  // Here hue_in_degrees is H, saturation is S_{HSV}, and value is V
  hue_in_degrees %= 360;

  // chroma is C; C = V \times S_{HSV}
  float chroma = saturation * value;
  // hue_sector is H'; H' = H / 60 degrees
  float hue_sector = hue_in_degrees / 60;
  float a;
  // a is H' mod 2
  for (a = hue_sector; a >= 2.0; a -= 2.0) 
    ;
  // a is H' mod 2 - 1
  a -= 1.0;
  // X = C(1 - |H' mod 2 - 1|)
  float X = chroma * (1 - (a >= 0.0 ? a : -a) );
  float red, green, blue;
  if (hue_sector < 1) {
    red = chroma;
    green = X;
    blue = 0.0;
  } else if (hue_sector < 2) {
    red = X;
    green = chroma;
    blue = 0.0;
  } else if (hue_sector < 3) {
    red = 0.0;
    green = chroma;
    blue = X;
  } else if (hue_sector < 4) {
    red = 0.0;
    green = X;
    blue = chroma;
  } else if (hue_sector < 5) {
    red = X;
    green = 0.0;
    blue = chroma;
  } else {
    red = chroma;
    green = 0.0;
    blue = X;
  }
  float m = value - chroma;
  red += m;
  green += m;
  blue += m;
  rgb.r = floor(red * 255 + 0.5);
  rgb.g = floor(green * 255 + 0.5);
  rgb.b = floor(blue * 255 + 0.5);
}
