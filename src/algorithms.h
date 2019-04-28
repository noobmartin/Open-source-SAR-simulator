#ifndef _algorithms_h
#define _algorithms_h

#include "common_types.h"

void gbp(matrix* pc_image, matrix* sar_image, radar_variables* variables);
void gbp_fft(matrix* sar_image, matrix* sar_image_fft, radar_variables* variables);
void radar_imager(matrix* scene, matrix* radar_image, matrix* chirp, radar_variables* variables);
void insert_waveform_in_scene(matrix* chirp, matrix* scene, radar_variables* variables);
void pulse_compress_signal(matrix* chirp, matrix* match, matrix* pc_waveform, radar_variables* variables);
void pulse_compress_image(matrix* radar_image, matrix* pc_image, matrix* match, radar_variables* variables);


#endif