#ifndef _algorithms_h
#define _algorithms_h

#include "common_types.h"

void gbp(matrix* data_matrix, radar_variables* variables);
void gbp_fft(matrix* data_matrix, radar_variables* variables);
void radar_imager(matrix* data_matrix, radar_variables* variables);
void insert_waveform_in_scene(matrix* data_matrix, radar_variables* variables);
void pulse_compress_signal(matrix* data_matrix, radar_variables* variables);
void pulse_compress_image(matrix* data_matrix, radar_variables* variables);


#endif