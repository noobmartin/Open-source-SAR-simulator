#include "sar_simulator.h"

void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output){
  int i;
  for(i = 0; i < kernel_length; i++){
    kernel[i] /= kernel_length;
  }
  fftw_plan fft = fftw_plan_dft_1d(kernel_length, kernel, output, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
  for(i = 0; i < kernel_length; i++){
    kernel[i] *= kernel_length;
  }
}

void chirp_generator(data_arrays* data, radar_variables* variables){
  double end_time = (double)variables->btproduct/variables->bandwidth;
  double chirp_rate = variables->bandwidth/(2*end_time);

  unsigned long int sample_frequency = 5*variables->bandwidth;
  unsigned long int time_steps = end_time*sample_frequency;
  variables->chirp_length = time_steps;
  variables->signal_distance = end_time*C;

  printf("Chirp signal distance: %lfm\n", variables->signal_distance);
  
  data->chirp_time_vector = malloc(time_steps*sizeof(double));
  data->chirp_signal = malloc(time_steps*sizeof(double complex));

  int i;
  double last_time = 0;
  for(i = 0; i < end_time*sample_frequency; i++){
    data->chirp_time_vector[i] = last_time;
    last_time += (double)1/sample_frequency;
  }

  for(i = 0; i < time_steps; i++){
    //data->chirp_signal[i] = cexp(-_Complex_I*2*PI*(variables->start_frequency*data->chirp_time_vector[i]+chirp_rate*data->chirp_time_vector[i]*data->chirp_time_vector[i]));
    data->chirp_signal[i] = cos( 2*PI*(variables->start_frequency*data->chirp_time_vector[i]+chirp_rate*data->chirp_time_vector[i]*data->chirp_time_vector[i] ))+_Complex_I*sin( 2*PI*(variables->start_frequency*data->chirp_time_vector[i]+chirp_rate*data->chirp_time_vector[i]*data->chirp_time_vector[i] ));
  }
}

void chirp_matched_generator(data_arrays* data, radar_variables* variables){
  data->matched_chirp = malloc(variables->chirp_length*sizeof(complex double));

  complex double* chirpfft = malloc(variables->chirp_length*sizeof(complex double));
  fftw_plan chirpfftp = fftw_plan_dft_1d(variables->chirp_length, data->chirp_signal, chirpfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(chirpfftp);
  int z;
  for(z = 0; z < variables->chirp_length; z++){
    chirpfft[z] = conj(chirpfft[z]);
  }

  fftw_plan matchedp = fftw_plan_dft_1d(variables->chirp_length, chirpfft, data->matched_chirp, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(matchedp);

  free(chirpfft);

  fftw_destroy_plan(chirpfftp);
  fftw_destroy_plan(matchedp);
}

float calculate_compressed_pulse_resolution(data_arrays* data, radar_variables* variables){
  double complex* waveform = data->pulse_compressed_waveform;

  // Find maximum amplitude position
  int i, max_position;
  double max_amplitude = 0;
  for(i = 0; i < variables->chirp_length; i++){
    if( sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2) ) > max_amplitude){
      max_amplitude = waveform[i];
      max_position = i;
    }
  }
  double half_amplitude = max_amplitude / 2;
  int low_index = 0;
  double low_value;
  for(i = max_position; i >= 0; i--){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	low_index = i;
	break;
    }
  }
  int high_index = 0;
  for(i = max_position; i < variables->chirp_length; i++){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	high_index = i;
	break;
    }
  }

  return variables->signal_distance*(high_index-low_index)/variables->chirp_length;
}
