/* Note to the reader:
 * To use the standard "complex" data type in C, complex.h must be included
 * before fftw3.h - otherwise a bunch of compiler warnings will arise.
 */
#include <complex.h>
#include <fftw3.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common_functions.h"
#include "waveforms.h"

void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output){
  for(int i = 0; i < kernel_length; i++){
    kernel[i] /= kernel_length;
  }

  fftw_plan fft = fftw_plan_dft_1d(kernel_length, kernel, output, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);

  for(int i = 0; i < kernel_length; i++){
    kernel[i] *= kernel_length;
  }

}

void chirp_generator(matrix* data, radar_variables* variables){
  printf("Chirp start frequency (Hz): ");
  int ret = scanf("%li", &variables->start_frequency);

  if(variables->start_frequency < 0){
    printf("Negative start frequency entered, closing.\n");
	  return;
  }

  printf("Chirp bandwidth (Hz): ");
  ret = scanf("%li", &variables->bandwidth);

  if(variables->bandwidth < 0){
	  printf("Negative bandwidth entered, closing.\n");
	  return;
  }

  printf("Chirp bandwidth-time product: ");
  ret = scanf("%u", &variables->btproduct);
    
  if(variables->btproduct < 1){
	  printf("Too small BT-product entered, closing.\n");
	  return;
  }

  double chirp_duration         = (double)variables->btproduct/variables->bandwidth;
  double chirp_rate             = variables->bandwidth/(2*chirp_duration);
  unsigned int sample_frequency = 2*variables->bandwidth;
  variables->chirp_length       = chirp_duration*sample_frequency;
  variables->signal_distance    = chirp_duration*C;

  printf("Chirp rate: %lf (Hz/s)\n",             chirp_rate);
  printf("Sample frequency: %u (Hz)\n",        sample_frequency);
  printf("Chirp samples: %lf\n",           variables->chirp_length);
  printf("Uncompressed pulse resolution: %fm\n", variables->signal_distance);
  
  /* Allocate memory for the time vector and generate it. */
  matrix* meta_time_vector = get_last_node(data);
  strcpy(meta_time_vector->name, "time_vector");
  meta_time_vector->rows = variables->chirp_length;
  meta_time_vector->cols = 1;
  meta_time_vector->data = malloc(variables->chirp_length*sizeof(complex double));

  double last_time = 0;
  for(int i = 0; i < chirp_duration*sample_frequency; i++){
    meta_time_vector->data[i] = last_time;
    last_time += (double)1/sample_frequency;
  }
  
  /* Allocate memory for the chirp waveform and generate it. */
  matrix* meta_chirp = get_last_node(data);
  strcpy(meta_chirp->name, "chirp");
  meta_chirp->rows = variables->chirp_length;
  meta_chirp->cols = 1;
  meta_chirp->data = malloc(variables->chirp_length*sizeof(complex double));

  for(int i = 0; i < variables->chirp_length; i++){
    double time = meta_time_vector->data[i];
    meta_chirp->data[i] = cexp(_Complex_I*2*PI*(variables->start_frequency*time+chirp_rate*time*time));
  }
  
}

void chirp_matched_generator(matrix* data, radar_variables* variables){
  /* Performs the following operations:
   * Chirp waveform -> FFT -> Complex conjucation -> Division by number of rows -> IFFT -> Matched chirp waveform
   */
  matrix* chirp    = get_matrix(data, "chirp");
  
  /* Execute FFT on the chirp waveform. */
  complex double* chirpfft  = malloc(chirp->rows*sizeof(complex double));
  fftw_plan       chirpfftp = fftw_plan_dft_1d(chirp->rows, chirp->data, chirpfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(chirpfftp);
 
  /* Complex-conjugate and divide by the number of rows. */
  for(int z = 0; z < chirp->rows; z++){
    chirpfft[z] = conj(chirpfft[z])/chirp->rows;
  }

  /* Allocate memory for the matched waveform */
  matrix* match = get_last_node(data);
  match->rows   = chirp->rows;
  match->cols   = chirp->cols;
  strcpy(match->name, "match");
  match->data = malloc(match->rows*sizeof(complex double));
  
  /* Execute IFFT to get the matched waveform. */
  fftw_plan matchedp = fftw_plan_dft_1d(chirp->rows, chirpfft, match->data, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(matchedp);

  free(chirpfft);

  fftw_destroy_plan(chirpfftp);
  fftw_destroy_plan(matchedp);
}

float calculate_compressed_pulse_resolution(matrix* data, radar_variables* variables){
  matrix* meta_pc_waveform = get_matrix(data, "pulse_compressed_waveform");
  double complex* waveform = meta_pc_waveform->data;

  // Find maximum amplitude position
  int    max_position  = 0;
  double max_amplitude = 0;
  int    low_index     = 0;
  int    high_index    = 0;

  for(int i = 0; i < meta_pc_waveform->rows; i++){
    if( sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2) ) > max_amplitude){
      max_amplitude = waveform[i];
      max_position = i;
    }
  }

  double half_amplitude = max_amplitude / 2;
  double low_value;
  for(int i = max_position; i >= 0; i--){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	    low_index = i;
	    break;
    }
  }

  for(int i = max_position; i < variables->chirp_length; i++){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	    high_index = i;
	    break;
    }
  }

  return variables->signal_distance*(high_index-low_index)/meta_pc_waveform->rows;
}
