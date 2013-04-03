/* 
gsl-sprng.h

Code declaring a new GSL random number type "gsl_rng_sprng20"
which is a thin wrapper over the SPRNG 2.0 parallel random
number generator.

To use, just add the line:
#include "gsl-sprng.h"
immediately after the line:
#include <gsl/gsl_rng.h>
near the start of your code. The new type should now be available.
Make sure you alloc the rng on each processor. If you wish to
set a seed, you should set it to be the same on each processor.

Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 28/8/2002

*/


#define SIMPLE_SPRNG
#define USE_MPI
#include "sprng.h"

static void sprng_set(void * vstate,unsigned long int s)
{
  init_sprng(DEFAULT_RNG_TYPE,s,SPRNG_DEFAULT);
}

static unsigned long sprng_get(void * vstate)
{
  return( (long) isprng() );
}

static double sprng_get_double(void * vstate)
{
  return( (double) sprng());
}

static const gsl_rng_type sprng_type =
  {"sprng20",        /* name */
   0x7fffffffUL,     /* RAND_MAX */
   0,                /* RAND_MIN */
   0,                /* size of state - not sure about this */
   &sprng_set,          /* initialisation */
   &sprng_get,          /* get integer RN */
   &sprng_get_double};  /* get double RN */

const gsl_rng_type *gsl_rng_sprng20 = &sprng_type;

/* eof */

