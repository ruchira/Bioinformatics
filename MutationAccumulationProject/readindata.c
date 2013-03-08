// readindata.c
// Authors: Trevor Graham, Ruchira S. Datta
// Copyright (c) 2012, Regents of the University of California
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

/* Reads data from input files and defines necessary structures to handle data */
/* Trevor Graham CoMPlEX, UCL, April 2005 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "readindata.h"

initialconditions ics;

////////////////////////////////////////////////////////////////////////////////////////////////
void readics( int argc, char* argv[] ) {
  // parses arguments from command line.  Uses processArg to extract data
  printf("#reading ICs... ");
      
  int i;
  unsigned int j;
  char param[ 100 ];      /* should be enough space for the param! */
  char value[ 200 ];      /* should be enough space for the value! */
  char found_equals;      /* Index of where the equals is */

  /* Store the program name for future ref. 
   * No need to copy argv[0] since it will always be in scope) */
  ics.progName = argv[ 0 ];
  /* Check if we just want to display the usage for the prog */
  if( argc == 1 )
      Usage( );

  /* set the default parameter values */
    
  strcpy(ics.filename, "out.txt");
  strcpy(ics.run_control_file, "");
  strcpy(ics.all_params, "");
  ics.run_number=0;
	ics.Ninitial = 1000;
	ics.repeats = 1;
	ics.gens = 1500;
  ics.seed = 3141592;
  ics.growthModel = 'K';
	ics.mostagg = 20;
  ics.muB = 0.001;
	ics.cancerProp = 0.1;
	ics.mutProp = 0.1;
  ics.Md = 50;
  ics.Mh = 50;
  ics.Mp = 50;
  ics.Mm = 50;
	ics.mut_max = 10;
  ics.num_passengers_needed_for_cancer = 0;
  ics.probD = 0.0001;
  ics.probH = 0.001;
  ics.probM = 0.00005;
  ics.sD = 0.01;
  ics.sH = 0.01;;
  ics.sM = 0.01;
  ics.beta = 0.01;
  ics.output_marginalized_counts = 0;
    
  // NOTE!!! This procedure does not do any range checking.  This is dangerous
  // and can cause a buffer overflow.

  /* Parse the program arguments one at a time.
   * Note - all args should be of the form -param=value */
  for( i = 1; i < argc; ++i )
  {
      /* Ensure we start with a '-' */
      if( argv[ i ][ 0 ] != '-' )
      {
          printf( "Argument %d should start with a '-'\n", i );
          exit(0);
      }

      /* Now find the equals */
      found_equals = 0;
      for( j = 1; j < strlen( argv[ i ] ); ++j )
      {
          if( argv[ i ][ j ] == '=' )
          {
              found_equals = j;
              param[ j - 1 ] = '\0';
          }
          else
          {
              /* If we haven't come across the '=' yet we are in the param
               * If we have already found the '=' we are in the value */
              if( found_equals == 0 )
              {
                  param[ j - 1 ] = argv[ i ][ j ];
              }
              else
              {
                  value[ j - found_equals - 1 ] = argv[ i ][ j ];
              }
          }
      }

      /* No '=' found therefore we have an invalid argument */
      if( found_equals == 0 )
      {
          printf( "Could not find '=' in arg %d\n", i );
          exit(0);
      }

      /* Ensure nul string termination then process this param, value pair */
      value[ j - found_equals - 1 ] = '\0';
      ProcessArg( param, value);
  }
  // If we do have a run control file, read in command-line arguments from
  // there.  These will override any that were on the invocation command line.

  if (ics.run_control_file[0] != '\0') {
    ReadInitialConditionsFromRunControlFile();
  }
  printf("done.\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessArg( const char* param, char* value)
{
    
    // parse value
    
    if( strcmp( param, "repeats" ) == 0 )
    {
        ics.repeats = atoi( value );
    }
    else if( strcmp( param, "gens" ) == 0 )
    {
        ics.gens = atoi( value );
    }
    else if( strcmp( param, "Ninitial" ) == 0 )
    {
        ics.Ninitial = atoi( value );
    }
    else if( strcmp( param, "seed" ) == 0 )
    {
        ics.seed = atoi( value );
    }
    else if( strcmp( param, "growthModel" ) == 0 )
    {
        ics.growthModel = *value;
    }
    else if( strcmp( param, "filename" ) == 0 )
    {
        strcpy(ics.filename, value);
    }
    else if( strcmp( param, "run_control_file" ) == 0 )
    {
        strcpy(ics.run_control_file, value);
    }
    else if( strcmp( param, "run_number" ) == 0 )
    {
      ics.run_number = atoi( value );
    }
    else if( strcmp( param, "beta" ) == 0 )
    {
        ics.beta = atof( value );
    }	
    else if( strcmp( param, "muB" ) == 0 )
    {
        ics.muB = atof( value );
    }
    else if( strcmp( param, "sH" ) == 0 )
    {
        ics.sH = atof( value );
    }
    else if( strcmp( param, "sD" ) == 0 )
    {
        ics.sD = atof( value );
    }
    else if( strcmp( param, "sM" ) == 0 )
    {
        ics.sM = atof( value );
    }
    else if( strcmp( param, "Md" ) == 0 )
    {
        ics.Md = atoi( value );
    }
    else if( strcmp( param, "Mh" ) == 0 )
    {
        ics.Mh = atoi( value );
    }
    else if( strcmp( param, "Mp" ) == 0 )
    {
        ics.Mp = atoi( value );
    }
    else if( strcmp( param, "Mm" ) == 0 )
    {
        ics.Mm = atoi( value );
    }
    else if( strcmp( param, "probD" ) == 0 )
    {
        ics.probD = atof( value );
    }
    else if( strcmp( param, "probH" ) == 0 )
    {
        ics.probH = atof( value );
    }
    else if( strcmp( param, "probM" ) == 0 )
    {
        ics.probM = atof( value );
    }
	else if( strcmp( param, "cancerProp") == 0)
	{
		ics.cancerProp = atof(value);
	}
	else if( strcmp( param, "mutProp") == 0)
	{
		ics.mutProp = atof(value);
	}
    else if( strcmp( param, "mostagg" ) == 0 )
    {
        ics.mostagg = atoi( value );
    }
    else if( strcmp( param, "mut_max" ) == 0 )
    {
        ics.mut_max = atoi( value );
    }	
    else if( strcmp( param, "num_passengers_needed_for_cancer" ) == 0 )
    {
        ics.num_passengers_needed_for_cancer = atoi( value );
    }	
    else if( strcmp( param, "output_marginalized_counts" ) == 0 )
    {
      ics.output_marginalized_counts = atoi( value );
    }
    else
    {
        printf( "Unrecognised program argument: <%s>\n", param );
        Usage( );
    }
    //printf( "Param = <%s>, Value = <%s>\n", param, value );
}


////////////////////////////////////////////////////////////////////////////////////////////////
void Usage( void )
{
    /* Displays the program usage and then exits */
    printf( "\n\nUsage:\n" );
    printf("\t-run_control_file=file from which to parse initial conditions\n");
    printf("\t-run_number=the task id of the current job\n");
    printf("\t-repeats=number repeats\n");
    printf("\t-gens=number generation per repeats\n");
    printf("\t-Ninitial=initial number of cells in tumour mass\n");
    printf("\tbeta=growth rate");
    printf("\tmuB=background mutation rate\n");
    printf("\tMd=number tumour suppresor/oncogenes\n");
    printf("\tMh=number housekeeper genes\n");
    printf("\tMp=number passenger loci\n");
    printf("\tMm=number mutator genes\n");
    printf("\tprobD=probability that a mutation is in a driver gene\n");
    printf("\tprobH=probability that a mutation is in a housekeeper gene\n");
    printf("\tprobM=probability that a mutation is in a mutator gene\n");
    printf("\tsD=selective advantage of a TS/oncogene\n");
    printf("\tsH=selective disadvantage of a housekeeper\n");
    printf("\tsM=increase in mutation rate of mutator gene\n");
    printf("\tmostagg=mutations to a cancer\n");
    printf("\tcancerProp=proportion of cancer that must have \"mostagg\" mutations for cancer to be called\n");
    printf("\tmutProp=proportion of cancer that must have a mutator phenotype for \"mutator_tumour\" to be called\n");
    //printf("\tmut_max=maximum number of simultaneous mutations in a single cell >=10\n");
    printf("\t-seed=random number generators seed\n");
    printf("\t-growthModel=population size as function of time\n");
    printf("\t            K (constant), Q (quadratic), C (cubic), E (exponential)\n");
    printf("\t-filename=output filename\n");    
    printf("NOTE: Parameter values in the run control file (if present) take precedence,\n");
    printf("then any command-line options given,\n");
    printf("and finally the default parameter values.\n");
    exit( 1 );
}

////////////////////////////////////////////////////////////////////////////////////////////////
void printParameters( void )
{
    /* prints the program parameters */
    printf("#***************************************\n");
    printf("#repeats=%d\n",ics.repeats);
    printf("#gens=%d\n",ics.gens);
    printf("#Ninitial=%d\n",ics.Ninitial);
    printf("#beta=%f\n",ics.beta);
    printf("#muB=%f\n",ics.muB);
    printf("#Md=%d\n",ics.Md);
    printf("#Mh=%d\n",ics.Mh);
    printf("#Mp=%d\n",ics.Mp);
    printf("#Mm=%d\n",ics.Mm);
    printf("#sD=%f\n",ics.sD);
    printf("#sH=%f\n",ics.sH);
    printf("#sM=%f\n",ics.sM);
    printf("#probD=%f\n",ics.probD);
    printf("#probH=%f\n",ics.probH);
    printf("#probM=%f\n",ics.probM);
    printf("#mostagg=%d\n",ics.mostagg);
    printf("#cancerProp=%f\n",ics.cancerProp);
    printf("#mutProp=%f\n",ics.mutProp);
    printf("#mut_max=%d\n",ics.mut_max);
    printf("#num_passengers_needed_for_cancer=%d\n",
            ics.num_passengers_needed_for_cancer);
    printf("#seed=%lu\n",ics.seed);
    printf("#growthModel=%c\n",ics.growthModel);
    printf("#output_marginalized_counts=%d\n", ics.output_marginalized_counts);
	printf("#filename=%s\n",ics.filename);
    printf("#run_control_file=%s\n", ics.run_control_file);
    printf("#all_params=%s\n", ics.all_params);
    printf("#run_number=%d\n", ics.run_number);
    printf("#***************************************\n");
}

void printParametersToFile(FILE *fp)
{
    /* prints the program parameters */
    fprintf(fp,"#***************************************\n");
    fprintf(fp,"#repeats=%d\n",ics.repeats);
    fprintf(fp,"#gens=%d\n",ics.gens);
	fprintf(fp,"#Ninitial=%d\n",ics.Ninitial);
    fprintf(fp,"#beta=%f\n",ics.beta);
    fprintf(fp,"#muB=%f\n",ics.muB);
    fprintf(fp,"#Md=%d\n",ics.Md);
    fprintf(fp,"#Mh=%d\n",ics.Mh);
    fprintf(fp,"#Mp=%d\n",ics.Mp);
    fprintf(fp,"#Mm=%d\n",ics.Mm);
    fprintf(fp,"#sD=%f\n",ics.sD);
    fprintf(fp,"#sH=%f\n",ics.sH);
    fprintf(fp,"#sM=%f\n",ics.sM);
    fprintf(fp,"#probD=%f\n",ics.probD);
    fprintf(fp,"#probH=%f\n",ics.probH);
    fprintf(fp,"#probM=%f\n",ics.probM);
    fprintf(fp,"#mostagg=%d\n",ics.mostagg);
    fprintf(fp,"#cancerProp=%f\n",ics.cancerProp);
    fprintf(fp,"#mutProp=%f\n",ics.mutProp);
    fprintf(fp,"#mut_max=%d\n",ics.mut_max);
    fprintf(fp,"#num_passengers_needed_for_cancer=%d\n",
            ics.num_passengers_needed_for_cancer);
    fprintf(fp,"#seed=%lu\n",ics.seed);
    fprintf(fp,"#growthModel=%c\n",ics.growthModel);
    fprintf(fp,"#output_marginalized_counts=%d\n", 
            ics.output_marginalized_counts);
	fprintf(fp,"#filename=%s\n",ics.filename);
    fprintf(fp,"#run_control_file=%s\n", ics.run_control_file);
    fprintf(fp,"#all_params=%s\n", ics.all_params);
    fprintf(fp,"#run_number=%d\n", ics.run_number);
    fprintf(fp,"#***************************************\n");
}

void ReadInitialConditionsFromRunControlFile( void ) {
  // The format of the run control file is a bunch of lines with sets of
  // command-line arguments
  // For instance:
  // -repeats 2000 -sD 0.1 -sH 0.01
  // -repeats 1000 -sD 0.01 -sH 0.1
  // -repeats 100 -sM 4e-5 -sD 0.01 -sH 0.01
  // This means: the first batch of tasks should do 2000 repeats each with sD
  // 0.1 and sH 0.01, the next should do 1000 repeats each with sD 0.01 and sH
  // 0.1, the should do 100 repeats each with sD 0.01, sH 0.01, and sM 4e-5.
  // This procedure goes through the run control file and finds the line
  // corresponding to the ics.run_number argument.  Then this procedure
  // processes that line just like command-line arguments.  The procedure also
  // modifies the output filename to reflect these specific parameter values,
  // and the run number.

  // NOTE!!! This procedure does not do any range checking on strings.  This is
  // dangerous and can cause a buffer overflow.
  FILE *fp;
  char buffer[256];
  char param[100];
  char value[200];
  int line_number, arg_number;
  char found_equals;      // Whether we have found the equals sign

  char *argptr, *cur_ptr, *param_ptr, *value_ptr, *all_params_ptr;


  if (!(fp = fopen(ics.run_control_file, "r"))) {
    printf("#*** error! couldn't open run control file %s for reading\n",
          ics.run_control_file);
    exit(0);
  }
  line_number = 0;
  while (fgets(buffer, 256, fp) != NULL) {
    ++line_number;
    if (line_number == ics.run_number) {
      // We have found the line number we want! Process this line.
      arg_number = 0;
      all_params_ptr = ics.all_params;
      // Start the parameter value specification string with a leading
      // underscore
      *all_params_ptr++ = '_';
      // Start parsing space-delimited arguments from the rest of the line
      for (argptr = strtok(buffer, " \n"); 
          argptr != NULL; 
          argptr = strtok(NULL, " \n")) {
        ++arg_number;
        cur_ptr = argptr;
        param_ptr = param;
        value_ptr = value;
        if( *cur_ptr != '-' )
        {
            printf( "Argument %d on line %d of %s should start with a '-'\n",
                      arg_number, line_number, ics.run_control_file );
            exit(0);
        }
        /* Now find the equals */
        found_equals = 0;
        // The initialization of this loop will move past the dash
        for (++cur_ptr; *cur_ptr; ++cur_ptr) {
          if (*cur_ptr == '=') {
            // We've come to the end of the parameter name
            found_equals = 1;
            // Terminate the parameter name
            *param_ptr = '\0';
            // Delimit the parameter name in the parameter value specification
            // string with an underscore
            *all_params_ptr++ = '_';
            // Get out of the for loop for reading the parameter name
            break;
          } else {
            // Keep copying the parameter name
            *param_ptr++ = *cur_ptr;
            *all_params_ptr++ = *cur_ptr;
          }
        }

        if (!found_equals || !*cur_ptr) {
          printf( "Could not find '=' in arg %d on line %d of %s\n", 
                  arg_number, line_number, ics.run_control_file );
          exit(0);
        }
        int i = 0;
        // The initialization of this loop will move past the equals sign
        for (++cur_ptr; *cur_ptr; ++cur_ptr) {
          ++i;
          // Keep copying the value string
          *value_ptr++ = *cur_ptr;
          *all_params_ptr++ = *cur_ptr;
        }
        // Terminate the value string
        *value_ptr = '\0';
        // Delimit the value in the parameter value specification string with
        // an underscore (note that this will leave a trailing underscore at
        // the end)
        *all_params_ptr++ = '_';

        // Process the argument
        ProcessArg( param, value );

      }
      // Terminate the parameter value specification string
      *all_params_ptr = '\0';
      // Loop to the end of ics.filename
      for (cur_ptr = ics.filename; *cur_ptr; ++cur_ptr) {
        ;
      }
      // Now cur_ptr is pointing to the '\0' ending ics.filename
      // Concatenate the parameter value specification and run number to the
      // end of the filename
      // NOTE!!! This is dangerous and can cause a buffer overflow
      sprintf(cur_ptr, "%srcline_%05d", ics.all_params, ics.run_number);
      // Stop reading the run control file
      break;
    }
  }
  // If we never find the current job number in the whole run control file, we
  // won't process any command line arguments from the file; we will just use
  // the initial conditions specified on this program's own command line, or
  // the default values.
  fclose(fp);
}
