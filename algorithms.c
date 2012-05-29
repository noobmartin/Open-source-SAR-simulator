/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 */

#include "sar_simulator.h"

/* GBP is the main SAR image formation (or de-skewness algorithm as I like to call it).
 * It focuses the hyperbolas which are a result of the stripmap radar scan into point-like objects with sidelobes.
 * In general, it does a superb job, if the input data is "nice enough".
 * It handles arbitrary scene sizes, and is thus only restricted by computer memory and processing time.
 * The downside of GBP is its processing time, which is in the order of O(N^3).
 * The normalization is necessary to avoid floating point exceptions and nan values later on in the processing stage.
 */
void gbp(data_arrays* data, radar_variables* variables){
  int i;
  for(i = 0; i < variables->ncols*variables->nrows; i++){
    if(isnan(data->pulse_compressed_radar_image[i]))
      data->pulse_compressed_radar_image[i] = 0;
  }

  data->sar_image = malloc(variables->nrows*variables->ncols*sizeof(double complex));
  int j,k,l;
  for(j = 0; j < variables->ncols; j++){
    for(k = 0; k < variables->nrows; k++){
      for(l = 0; l < variables->ncols; l++){
	unsigned int range_index = sqrt((l-j)*(l-j)+k*k);
	if(range_index < variables->nrows){
	  data->sar_image[j*variables->nrows+k] += data->pulse_compressed_radar_image[l*variables->nrows+range_index]/variables->nrows;
	  if(isnan(data->sar_image[j*variables->nrows+k])){
	    data->sar_image[j*variables->nrows+k] = 0;
	  }
	}
      }
    }
  }

  normalize_image(data->sar_image, variables->nrows, variables->ncols);

}

/* To evaluate the quality of a simulated SAR image, one often looks at the spectrum of the GBP image.
 * For a point object in the GBP image, the 2D FFT should have a "fan" shape.
 */
void gbp_fft(data_arrays* data, radar_variables* variables){
  data->sar_fft = malloc(variables->nrows*variables->ncols*sizeof(double complex));
  data->sar_img_shifted = malloc(variables->nrows*variables->ncols*sizeof(double complex));
  // This shifts the frequencies of the 2D FFT so that zero frequency is in the middle of the image.
  int i, j;
  for(i = 0; i < variables->ncols; i++){
    for(j = 0; j < variables->nrows; j++){
	data->sar_img_shifted[variables->nrows*i+j] = data->sar_image[variables->nrows*i+j]*pow(-1,i+j);
    }
  }

  fftw_plan fft = fftw_plan_dft_2d(variables->ncols, variables->nrows, data->sar_img_shifted, data->sar_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
}

/* This algorithm is designed by me, and simulates scanning over a scene.
 * It takes antenna azimuth beamwidth and antenna elevation into account, and calculates how many azimuth bins the antenna sees from each azimuth position.
 * All the data within the antenna's field of view is compressed into a single column, for each column, taking distance into account.
 * This is what generates the familiar hyperbola in the raw radar image.
 */
void radar_imager(data_arrays* data, radar_variables* variables){
  if(!data->scene_with_waveform){
    data->scene_with_waveform = malloc(variables->ncols*variables->nrows*sizeof(complex double));
    memset(data->scene_with_waveform, 0, variables->ncols*variables->nrows*sizeof(complex double));
  }
  
  data->undistorted_radar_image = malloc(variables->nrows*variables->ncols*sizeof(double complex));
  double azimuth_coverage = round(variables->altitude*tan(0.5*variables->beamwidth));
  unsigned int beamcrossrange = variables->chirp_length*azimuth_coverage/variables->signal_distance;

  printf("Antenna beam covers %f meters, %i columns.\n", azimuth_coverage, beamcrossrange);

  unsigned int i,j,k;
  unsigned int dist;
  int beam_value;
  for(i = 0; i < variables->ncols; i++){
    for(j = 0; j < 2*beamcrossrange; j++){
      beam_value = j-beamcrossrange;
      for(k = 0; k < variables->nrows; k++){
	if(i+beam_value < variables->ncols){
	  unsigned int dist = sqrt(pow(beam_value,2)+pow(k, 2));
	  if(dist < variables->nrows){
	    if(i+beam_value >= 0){
	      //data->undistorted_radar_image[i*variables->nrows+dist] += (data->scene_with_waveform[(i+beam_value)*variables->nrows+k]/pow(dist,2));
	      data->undistorted_radar_image[i*variables->nrows+dist] += (data->scene_with_waveform[(i+beam_value)*variables->nrows+k]);
	    }
	  }
	}
      }
    }
  }

  if(!data->radar_image){
    data->radar_image = malloc(variables->nrows*variables->ncols*sizeof(complex double));
    memcpy(data->radar_image, data->undistorted_radar_image, variables->nrows*variables->ncols*sizeof(complex double));
  }
}

// This just copies the waveform into the center of the scene, to simulate an object at the center.
void insert_waveform_in_scene(data_arrays* data, radar_variables* variables){
	data->scene_with_waveform = malloc(variables->ncols*variables->nrows*sizeof(double complex));
	int i;
	for(i = 0; i < variables->chirp_length; i++){
	  data->scene_with_waveform[ (int)(variables->ncols/2)*variables->nrows + (int)(variables->nrows/2) -variables->chirp_length/2 +i] = data->chirp_signal[i];
	  //data->scene_with_waveform[ (int)(variables->ncols/4)*variables->nrows + (int)(variables->nrows/2) -variables->chirp_length/2 +i] = data->chirp_signal[i];
	}
}

/* Pulse compression of a single waveform (it's nice to look at the plots to see that pulse compression works).
 * This is just convolution computed in the frequency domain with the FFT.
 */
void pulse_compress_signal(data_arrays* data, radar_variables* variables){
  data->pulse_compressed_waveform = malloc(variables->chirp_length*sizeof(double complex));
  int kernel_length = variables->chirp_length;
  unsigned int filter_length = 2*kernel_length;

  fftw_complex* padded_signal = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  data->pulse_compressed_waveform = malloc(kernel_length*sizeof(double complex));
  
  fftw_complex* sigfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* matchfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));
  
  fftw_plan sigp = fftw_plan_dft_1d(filter_length, padded_signal, sigfft, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan matchp = fftw_plan_dft_1d(filter_length, padded_kernel, matchfft, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan iff = fftw_plan_dft_1d(kernel_length, product, data->pulse_compressed_waveform, FFTW_FORWARD, FFTW_MEASURE);
  
  memset(padded_signal, 0, filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_signal, data->chirp_signal, kernel_length*sizeof(fftw_complex));
  memcpy(padded_kernel, data->matched_chirp, kernel_length*sizeof(fftw_complex));

  fftw_execute(sigp);
  fftw_execute(matchp);

  int i;
  for(i = 0; i < filter_length; i++)
    product[i] = sigfft[i]*matchfft[i];

  fftw_execute(iff);

  fftw_destroy_plan(sigp);
  fftw_destroy_plan(matchp);
  fftw_destroy_plan(iff);

  fftw_free(sigfft);
  fftw_free(matchfft);
  fftw_free(product);
}

/* Convolution of the entire radar image with the matched waveform, again implemented with the FFT.
 */
void pulse_compress_image(data_arrays* data, radar_variables* variables){
  data->pulse_compressed_radar_image = malloc(variables->nrows*variables->ncols*sizeof(double complex));
  
  // Make sure input has valid values.
  int kernel_length = variables->chirp_length;
  int z;
  for(z = 0; z < kernel_length; z++){
    if(isnan(data->matched_chirp[z]))
      data->matched_chirp[z] = 0;
  }

  normalize_image(data->radar_image, variables->nrows, variables->ncols);

  unsigned int filter_length = variables->nrows + kernel_length;

  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* kernel_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memset(kernel_fft, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_kernel, data->matched_chirp, kernel_length*sizeof(fftw_complex));

  // Compute fft of filter kernel.
  fftw_plan kernelfft = fftw_plan_dft_1d(filter_length, padded_kernel, kernel_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(kernelfft);
  fftw_complex* output_column;
  fftw_complex* padded_column = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_column_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  int i,j,k;
  for(i = 0; i < variables->ncols; i++){
    double complex* column = &data->radar_image[i*variables->nrows];

    // Make sure we have valid values.
    for(k = 0; k < variables->nrows; k++){
	if(isnan(column[k]))
	  column[k] = 0;
    }

    output_column = &data->pulse_compressed_radar_image[i*variables->nrows];
    memset(padded_column, 0, filter_length*sizeof(fftw_complex));
    memcpy(padded_column, column, variables->nrows*sizeof(fftw_complex));
    memset(padded_column_fft, 0, filter_length*sizeof(fftw_complex));
    fftw_plan colfft = fftw_plan_dft_1d(filter_length, padded_column, padded_column_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(colfft);


    for(j = 0; j < filter_length; j++){
	product[j] = padded_column_fft[j]*kernel_fft[j];
    }

    fftw_plan colifft = fftw_plan_dft_1d(variables->nrows, product, output_column, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(colifft);

    fftw_destroy_plan(colfft);
    fftw_destroy_plan(colifft);

    // Re-normalize ifft output.
    normalize_image(output_column, variables->nrows, 1);
  }

  fftw_free(product);
  fftw_free(padded_column_fft);
  fftw_free(padded_column);
  fftw_free(padded_kernel);

  fftw_destroy_plan(kernelfft);
}
