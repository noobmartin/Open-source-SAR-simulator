/* Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This header file contains basic structures for managing all data structures in the simulator.
 */

#ifndef sar_simulator_h
#define sar_simulator_h

#define PI 3.14159265
#define C 300000000
#define FILE_OUTPUT_PRECISION "%.10f\t"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <fftw3.h>
#include <stdbool.h>

typedef struct{
  bool         is_complex;
  unsigned int rows;
  unsigned int cols;
}radar_metadata;

typedef enum{
  Simulate,
  Process
}mode_type;

// Just as the data_arrays structure, this struct keeps variables during the simulation run which is passed to member functions.
typedef struct{
  long unsigned int start_frequency;
  long unsigned int bandwidth;
  unsigned int      chirp_length;
  unsigned int      btproduct;
  int               altitude;
  float             beamwidth;
  double            signal_distance;
  mode_type         mode;
  char              radar_data_filename[255];
}radar_variables;

typedef struct matrix{
  complex double* data;
  unsigned int    rows;
  unsigned int    cols;
  char            name[255];
  struct matrix*  next;
}matrix;

/* Function prototypes.
 * Although it may look messy, every function has the variables and data structure pointers as argument.
 * It keeps data encapsulated since I don't like to declare global variables.
 * This could be made more "decent" if the simulator was to be re-written in c++.
 */
void chirp_generator(matrix* data_matrix, radar_variables* variables);
void chirp_matched_generator(matrix* data_matrix, radar_variables* variables);
void insert_waveform_in_scene(matrix* data_matrix, radar_variables* variables);
void radar_imager(matrix* data_matrix, radar_variables* variables);
void pulse_compress_image(matrix* data_matrix, radar_variables* variables);
void gbp(matrix* data_matrix, radar_variables* variables);
void gbp_fft(matrix* data_matrix, radar_variables* variables);
void pulse_compress_signal(matrix* data_matrix, radar_variables* variables);
void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output);
void filter_dc(matrix* data_matrix, radar_variables* variables);
int write_data(matrix* data_matrix, radar_variables* variables);
int simulate(matrix* data_matrix, radar_variables* variables);
void process_data(matrix* data_matrix, radar_variables* variables);
int read_radar_file(matrix* data_matrix, radar_variables* variables);
float calculate_compressed_pulse_resolution(matrix* data_matrix, radar_variables* variables);
void free_memory(matrix* data_matrix);
void normalize_image(complex double* image, unsigned int rows, unsigned int cols);
matrix* get_matrix(matrix* data, const char* name);
matrix* get_last_node(matrix* data);
void build_metadata(matrix* data, radar_variables* variables);

#endif
