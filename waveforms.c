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

  double end_time = (double)variables->btproduct/variables->bandwidth;
  double chirp_rate = variables->bandwidth/(2*end_time);

  unsigned int sample_frequency = 80*variables->bandwidth;

  unsigned long int time_steps = end_time*sample_frequency;
  variables->chirp_length = time_steps;
  variables->signal_distance = end_time*C;

  printf("Chirp signal distance: %lfm\n", variables->signal_distance);
 
  matrix* meta_time_vector = get_last_node(data);
  strcpy(meta_time_vector->name, "time_vector");
  meta_time_vector->rows = time_steps;
  meta_time_vector->cols = 1;
  complex double* chirp_time_vector = malloc(time_steps*sizeof(complex double));
  meta_time_vector->data = chirp_time_vector;

  matrix* meta_chirp = get_last_node(data);
  strcpy(meta_chirp->name, "chirp");
  meta_chirp->rows = time_steps;
  meta_chirp->cols = 1;
  complex double* chirp = malloc(time_steps*sizeof(complex double));
  meta_chirp->data = chirp;

  int i;
  double last_time = 0;
  for(i = 0; i < end_time*sample_frequency; i++){
    chirp_time_vector[i] = last_time;
    last_time += (double)1/sample_frequency;
  }

  for(i = 0; i < time_steps; i++){
    double time = chirp_time_vector[i];
    chirp[i] = cexp(_Complex_I*2*PI*(variables->start_frequency*time+chirp_rate*time*time));
  }
  
  printf("Uncompressed pulse resolution: %fm\n", variables->signal_distance);
}

void chirp_matched_generator(matrix* data, radar_variables* variables){
  matrix* meta_chirp = get_matrix(data, "chirp");
  complex double* chirp = meta_chirp->data;
  
  matrix* meta_match = get_last_node(data);
  meta_match->rows = meta_chirp->rows;
  meta_match->cols = meta_chirp->cols;
  strcpy(meta_match->name, "match");
  complex double* match = malloc(meta_match->rows*sizeof(complex double));
  meta_match->data = match;

  complex double* chirpfft = malloc(meta_chirp->rows*sizeof(complex double));
  fftw_plan chirpfftp = fftw_plan_dft_1d(meta_chirp->rows, chirp, chirpfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(chirpfftp);
  int z;
  for(z = 0; z < meta_chirp->rows; z++){
    chirpfft[z] = conj(chirpfft[z])/meta_chirp->rows;
  }

  fftw_plan matchedp = fftw_plan_dft_1d(meta_chirp->rows, chirpfft, match, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(matchedp);

  free(chirpfft);

  fftw_destroy_plan(chirpfftp);
  fftw_destroy_plan(matchedp);
}

float calculate_compressed_pulse_resolution(matrix* data, radar_variables* variables){
  matrix* meta_pc_waveform = get_matrix(data, "pulse_compressed_waveform");
  double complex* waveform = meta_pc_waveform->data;

  // Find maximum amplitude position
  int i, max_position = 0;
  double max_amplitude = 0;
  for(i = 0; i < meta_pc_waveform->rows; i++){
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

  return variables->signal_distance*(high_index-low_index)/meta_pc_waveform->rows;
}
