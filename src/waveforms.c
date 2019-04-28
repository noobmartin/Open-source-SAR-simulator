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

void chirp_generator(matrix* time_vector, matrix* chirp, radar_variables* variables){
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

  double chirp_duration_s          = (double)variables->btproduct/variables->bandwidth;
  double chirp_rate_hz_per_s       = variables->bandwidth/(chirp_duration_s);
  unsigned int sample_frequency_hz = 2*variables->bandwidth;
  double sample_period             = 1/(double)sample_frequency_hz;
  variables->chirp_samples         = chirp_duration_s*sample_frequency_hz;
  variables->signal_distance       = chirp_duration_s*C;

  printf("Chirp duration: %f (s)\n",             chirp_duration_s);
  printf("Chirp rate: %lf (Hz/s)\n",             chirp_rate_hz_per_s);
  printf("Sample frequency: %u (Hz)\n",          sample_frequency_hz);
  printf("Chirp samples: %lf\n",                 variables->chirp_samples);
  printf("Uncompressed pulse resolution: %fm\n", variables->signal_distance);
  
  /* Allocate memory for the time vector and generate it. */
  time_vector->rows = variables->chirp_samples;
  time_vector->cols = 1;
  time_vector->data = malloc(variables->chirp_samples*sizeof(complex double));

  double last_time = 0;
  for(int i = 0; i < chirp_duration_s*sample_frequency_hz; i++){
    time_vector->data[i] = last_time;
    last_time += sample_period;
  }
  
  /* Allocate memory for the chirp waveform and generate it. */
  chirp->rows = variables->chirp_samples;
  chirp->cols = 1;
  chirp->data = malloc(variables->chirp_samples*sizeof(complex double));

  for(int i = 0; i < variables->chirp_samples; i++){
    double time = time_vector->data[i];
    chirp->data[i] = cexp(_Complex_I*2*PI*(variables->start_frequency*time+chirp_rate_hz_per_s*time*time));
  }
  
}

void chirp_matched_generator(matrix* chirp, matrix* match){
  /* Performs the following operations:
   * Chirp waveform -> FFT -> Complex conjucation -> Division by number of rows -> IFFT -> Matched chirp waveform
   */
  
  /* Execute FFT on the chirp waveform. */
  complex double* chirpfft  = malloc(chirp->rows*sizeof(complex double));
  fftw_plan       chirpfftp = fftw_plan_dft_1d(chirp->rows, chirp->data, chirpfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(chirpfftp);
 
  /* Complex-conjugate and divide by the number of rows. */
  for(int z = 0; z < chirp->rows; z++){
    chirpfft[z] = conj(chirpfft[z])/chirp->rows;
  }

  /* Allocate memory for the matched waveform */
  match->rows   = chirp->rows;
  match->cols   = chirp->cols;
  match->data = malloc(match->rows*sizeof(complex double));
  
  /* Execute IFFT to get the matched waveform. */
  fftw_plan matchedp = fftw_plan_dft_1d(chirp->rows, chirpfft, match->data, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(matchedp);

  free(chirpfft);

  fftw_destroy_plan(chirpfftp);
  fftw_destroy_plan(matchedp);
}

float calculate_compressed_pulse_resolution(matrix* pc_waveform, radar_variables* variables){
  double complex* waveform = pc_waveform->data;

  // Find maximum amplitude position
  int    max_position  = 0;
  double max_amplitude = 0;
  int    low_index     = 0;
  int    high_index    = 0;

  for(int i = 0; i < pc_waveform->rows; i++){
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

  for(int i = max_position; i < variables->chirp_samples; i++){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	    high_index = i;
	    break;
    }
  }

  return variables->signal_distance*(high_index-low_index)/pc_waveform->rows;
}
