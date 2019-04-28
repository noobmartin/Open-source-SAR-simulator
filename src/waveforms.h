/*
 * Author:  Alexander Rajula
 * Contact: alexander@rajula.se
 */
 
#ifndef _waveforms_h
#define _waveforms_h
 
#include "common_types.h"
 
void chirp_generator(matrix* data_matrix, radar_variables* variables);
void chirp_matched_generator(matrix* data_matrix, radar_variables* variables);
float calculate_compressed_pulse_resolution(matrix* data_matrix, radar_variables* variables);
void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output);

#endif