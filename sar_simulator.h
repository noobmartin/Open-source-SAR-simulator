/* Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This header file contains basic structures for managing all data structures in the simulator.
 */

#ifndef sar_simulator_h
#define sar_simulator_h

#define PI 3.14159265
#define C 300000000
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <fftw3.h>

// This struct is passed around in the simulator for ease of access to all simulator functions, as well as keeping the data encapsulated.
typedef struct{
  double* chirp_time_vector;
  double* matched_time_vector;
  double complex* chirp_fft;
  double complex* matched_fft;
  double complex* pulse_compressed_waveform;
  double complex* chirp_signal;
  double complex* matched_chirp;
  double complex* scene_with_waveform;
  double complex* undistorted_radar_image;
  double complex* radar_image;
  double complex* unfiltered_radar_image;
  double complex* pulse_compressed_radar_image;
  double complex* sar_image;
  double complex* sar_fft;
  double complex* sar_img_shifted;
  double complex* k_filter;
  double complex* phi_filter;
}data_arrays;

typedef struct{
  char real_or_complex;
  unsigned int rows;
  unsigned int cols;
}radar_metadata;

// Just as the data_arrays structure, this struct keeps variables during the simulation run which is passed to member functions.
typedef struct{
  long unsigned int start_frequency;
  long unsigned int bandwidth;
  unsigned int chirp_length;
  unsigned int nrows;
  unsigned int ncols;
  unsigned int btproduct;
  int altitude;
  float beamwidth;
  double signal_distance;
  char mode;
  char radar_data_filename[255];
  char real_or_complex_simulation[2];
  unsigned int k_filter_length;
  unsigned int phi_filter_length;
}radar_variables;

// This is only used for the apodization step.
typedef enum{
  HANNING,
  HAMMING,
  BARTLETT
}FILTER_TYPE;

/* Function prototypes.
 * Although it may look messy, every function has the variables and data structure pointers as argument.
 * It keeps data encapsulated since I don't like to declare global variables.
 * This could be made more "decent" if the simulator was to be re-written in c++.
 */
void chirp_generator(data_arrays* data, radar_variables* variables);
void chirp_matched_generator(data_arrays* data, radar_variables* variables);
void insert_waveform_in_scene(data_arrays* data, radar_variables* variables);
void radar_imager(data_arrays* data, radar_variables* variables);
void pulse_compress_image(data_arrays* data, radar_variables* variables);
void gbp(data_arrays* data, radar_variables* variables);
void gbp_fft(data_arrays* data, radar_variables* variables);
void pulse_compress_signal(data_arrays* data, radar_variables* variables);
void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output);
void filter_dc(data_arrays* data, radar_variables* variables);
int write_complex_data(data_arrays* data, radar_variables* variables);
int write_real_data(data_arrays* data, radar_variables* variables);
int simulate(data_arrays* data, radar_variables* variables);
void process_data(data_arrays* data, radar_variables* variables);
int read_radar_file(data_arrays* data, radar_variables* variables);
float calculate_compressed_pulse_resolution(data_arrays* data, radar_variables* variables);
void free_memory(data_arrays* data);
void generate_k_filter(data_arrays* data, radar_variables* variables, FILTER_TYPE type);
void generate_phi_filter(data_arrays* data, radar_variables* variables, FILTER_TYPE type);
void apodize_sar_fft(data_arrays* data, radar_variables* variables);
void rfi_suppression(data_arrays* data, radar_variables* variables);
void generate_gsm_interference(data_arrays* data, radar_variables* variables);
void normalize_image(complex double* image, unsigned int rows, unsigned int cols);

#endif
