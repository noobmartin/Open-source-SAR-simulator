#ifndef _common_types_h
#define _common_types_h

#include <stdbool.h>
#include <complex.h>

#define PI 3.14159265
#define C 300000000

typedef struct{
  bool         is_complex;
  unsigned int rows;
  unsigned int cols;
}radar_metadata;

// Just as the data_arrays structure, this struct keeps variables during the simulation run which is passed to member functions.
typedef struct{
  long unsigned int start_frequency;
  long unsigned int bandwidth;
  unsigned long int chirp_samples;
  unsigned int      btproduct;
  int               altitude;
  float             beamwidth;
  double            signal_distance;
  char              radar_data_filename[255];
}radar_variables;

typedef struct matrix{
  complex double* data;
  unsigned int    rows;
  unsigned int    cols;
  char            name[255];
}matrix;

#endif