/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 */

#include "sar_simulator.h"

/* This employs an adaptive line enhancer to the radar image.
 * Its purpose is to filter out radio frequency interference, for instance GSM and television broadcasts.
 * Tha ALE is controlled by the beta, epsilon, n_naught, and filter length parameters, and may be adjusted to achieve the desired effect.
 * Other RFI filtering approaches may be implemented and tested.
 */
/*
void rfi_suppression(data_arrays* data, radar_variables* variables){
  double beta = 0.3;
  double epsilon = 100000000;
  int n_naught = 1;
  unsigned int filter_length = 1024;
  double w[filter_length];
  memset(w, 0, filter_length*sizeof(double));

  //normalize_image(data->radar_image, variables->nrows, variables->ncols);

  double complex* filtered_image;
  filtered_image = malloc(variables->ncols*variables->nrows*sizeof(complex double));
  if(filtered_image == 0){
    return;
  }

  memset(filtered_image, 0, variables->ncols*variables->nrows*sizeof(complex double));

  int i,j;
  for(i = 0; i < variables->ncols*variables->nrows; i++){
	// Convolution
	for(j = 0; j < filter_length; j++){
	  if(i-n_naught-j >= 0)
	    filtered_image[i] += w[j]*(data->radar_image[i-n_naught-j]); 
	  if(isnan(filtered_image[i]))
	    filtered_image[i] = 0;
	}

	// Update filter coefficients
	for(j = 0; j < filter_length; j++){
	  if(j > 0 )
	    if(i >= n_naught)
	      w[j] = w[j-1] + ( (beta*conj(data->radar_image[i-n_naught])*(data->radar_image[i]-filtered_image[i]))/(pow(abs(data->radar_image[i-n_naught]),2)+epsilon));
	    else
	      w[j] = 0;
	  else
	    w[j] = 0;
	}
  }

  //normalize_image(filtered_image, variables->nrows, variables->ncols);

  // Save un-filtered image.
  data->unfiltered_radar_image = malloc(variables->ncols*variables->nrows*sizeof(complex double));
  memcpy(data->unfiltered_radar_image, data->radar_image, variables->ncols*variables->nrows*sizeof(complex double));

  memcpy(data->radar_image, filtered_image, variables->ncols*variables->nrows*sizeof(complex double));

  free(filtered_image);
}
*/

// Not fully implemented yet - ignore this function and do not run it during the simulation.
/*
void apodize_sar_fft(data_arrays* data, radar_variables* variables){
  generate_phi_filter(data, variables, HANNING);
  generate_k_filter(data, variables, HANNING);

  // Allocate memory for non-separable filter.
  complex double* window = malloc(variables->ncols*variables->nrows*sizeof(complex double));
  memset(window, 0, variables->ncols*variables->nrows*sizeof(complex double));
  
  double phi_naught = 2*atan(variables->ncols/(2*variables->nrows));
  unsigned int k_max = variables->nrows;
  unsigned int k_min = variables->nrows/2;
  unsigned int k_c = k_min+(k_max-k_min)/2;

  // Create filter from phi and k filters.
  int i,j;
  for(i = 0; i < variables->ncols; i++){
    for(j = 0; j < variables->nrows; j++){
      // This is not correct. It will crash.
      unsigned int atan_arg = ( i-(variables->ncols/2) ) / ( j-(variables->nrows/2) );
      unsigned int sqrt_arg = pow(i-(variables->ncols/2),2) + pow(j-(variables->nrows/2),2);
      unsigned int phi_arg;
      unsigned int k_arg;
      if(isinf(atan_arg)){
        phi_arg = 2*PI/phi_naught;
	k_arg = sqrt(sqrt_arg-k_c/(k_max-k_min));
	if(phi_arg > variables->phi_filter_length)
	  window[i*variables->nrows +j] = 0;
	else if(k_arg > variables->k_filter_length)
	  window[i*variables->nrows +j] = 0;
	else
	  window[i*variables->nrows+j] = data->phi_filter[phi_arg]*data->k_filter[k_arg];
      }
      else{
        phi_arg = 2*abs(atan(atan_arg))/phi_naught;
	k_arg = sqrt(sqrt_arg-k_c/(k_max-k_min));
	printf("phi_arg: %i\n", phi_arg);
	printf("k_arg: %i\n", k_arg);
	if(phi_arg > variables->phi_filter_length)
	  window[i*variables->nrows +j] = 0;
	else if(k_arg > variables->k_filter_length)
	  window[i*variables->nrows +j] = 0;
	else
          window[i*variables->nrows+j] = data->phi_filter[phi_arg]*data->k_filter[k_arg];
      }
    }
  }

  // Allocate memory for 2D FFT of filter.
  fftw_complex* window_fft = fftw_malloc(variables->ncols*variables->nrows*sizeof(fftw_complex));
  memset(window_fft, 0, variables->ncols*variables->nrows*sizeof(fftw_complex));

  // Calculate the 2D FFT of window.
  fftw_plan window_fft_plan = fftw_plan_dft_2d(variables->ncols, variables->nrows, window, window_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(window_fft_plan);
  fftw_destroy_plan(window_fft_plan);

  // Allocate memory for 2D FFT of 2D FFT of SAR image.
  fftw_complex* sar_fft_fft = fftw_malloc(variables->ncols*variables->nrows*sizeof(fftw_complex));
  memset(sar_fft_fft, 0, variables->ncols*variables->nrows*sizeof(fftw_complex));

  // Calculate 2D FFT of 2D FFT of SAR image.
  fftw_plan sar_fft_fft_plan = fftw_plan_dft_2d(variables->ncols, variables->nrows, data->sar_fft, sar_fft_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(sar_fft_fft_plan);
  fftw_destroy_plan(sar_fft_fft_plan);

  // CONVOLVE!
  for(i = 0; i < variables->ncols; i++){
    for(j = 0; j < variables->nrows; j++){
      sar_fft_fft[i*variables->nrows+j] *= window_fft[i*variables->nrows+j];
    }
  }

  // Inverse transform!
  fftw_plan sar_fft_ifft_plan = fftw_plan_dft_2d(variables->ncols, variables->nrows, sar_fft_fft, data->sar_fft, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(sar_fft_ifft_plan);
  fftw_destroy_plan(sar_fft_ifft_plan);

  // Don't stop there, get back the image as well!
  fftw_plan sar_ifft_plan = fftw_plan_dft_2d(variables->ncols, variables->nrows, data->sar_fft, data->sar_image, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(sar_ifft_plan);
  fftw_destroy_plan(sar_ifft_plan);
  
}
*/

// Generate filter coefficients for the apodization.
/*
void generate_k_filter(data_arrays* data, radar_variables* variables, FILTER_TYPE type){
  variables->k_filter_length = variables->nrows;
  data->k_filter = malloc(variables->k_filter_length*sizeof(complex double));
  
  int i;
  switch(type){
    case HANNING:
  	for(i = 0; i < variables->k_filter_length; i++){
    	  data->k_filter[i] = 0.5*(1-cos(2*PI*i/(variables->k_filter_length -1)));
  	}
    	break;
    case HAMMING:
	for(i = 0; i < variables->k_filter_length; i++){
	  data->k_filter[i] = 0.54 - 0.46*cos(2*PI*i/(variables->k_filter_length -1));
	}
    	break;
    case BARTLETT:
	for(i = 0; i < variables->k_filter_length; i++){
	  data->k_filter[i] = (2/(variables->k_filter_length -1))*((variables->k_filter_length-1)/2- abs(i-(variables->k_filter_length-1)/2));
	}
    	break;
  }
}
*/

// Generate filter coefficients for the apodization.
/*
void generate_phi_filter(data_arrays* data, radar_variables* variables, FILTER_TYPE type){
  variables->phi_filter_length = variables->nrows;
  data->phi_filter = malloc(variables->phi_filter_length*sizeof(complex double));

  int i;
  switch(type){
    case HANNING:
  	for(i = 0; i < variables->phi_filter_length; i++){
    	  data->phi_filter[i] = 0.5*(1-cos(2*PI*i/(variables->phi_filter_length -1)));
  	}
  	break;
    case HAMMING:
	for(i = 0; i < variables->phi_filter_length; i++){
	  data->phi_filter[i] = 0.54 - 0.46*cos(2*PI*i/(variables->phi_filter_length -1));
	}
    	break;
    case BARTLETT:
	for(i = 0; i < variables->phi_filter_length; i++){
	  data->phi_filter[i] = (2/(variables->phi_filter_length -1))*((variables->phi_filter_length-1)/2- abs(i-(variables->phi_filter_length-1)/2));
	}
    	break;
  }
}
*/

/* Although the GSM signal is continuous in nature, this "filter" mimics a stochastic GSM signal.
 * The generated GSM signal is simply added to the radar image.
 */
/*
void generate_gsm_interference(data_arrays* data, radar_variables* variables){
	// b[n] -> [P] -> r[n] -> [G] -> g[n] -> [S] -> phi[n] -> [cos] -> y[n]

	complex double* bsig = malloc(variables->nrows*variables->ncols*sizeof(complex double));
	if(!bsig)
	  return;

	srand(time(NULL));

	// Generate random binary signal.
	int i;
	for(i = 0; i < variables->nrows*variables->ncols; i++){
	  bsig[i] = rand()%2;
	}

	// Process pulse-train by Gaussian filter
	//double a_naught = 1;
	//double beta = 1;
	double f_c = 900000000;
	double t_naught = 0;
	double f_m = 250000;
	double tb = 0.3; // For GSM the time-bandwidth product is 0.3 for the Gaussian filter.
	double b = 0.387962624;
	double a = b/tb;
	unsigned int gaussian_filter_length = 1024;

	double* gaussian_filter = malloc(gaussian_filter_length*sizeof(double));
	if(!gaussian_filter)
	  return;

	for(i = 0; i < gaussian_filter_length; i++){
	  //gaussian_filter[i] = (a_naught/sqrt(PI))*beta*exp(-pow(i*beta, 2));
	  gaussian_filter[i] = (sqrt(PI)/a)*exp(-pow(PI*i/a,2));
	}

	unsigned int filter_length = gaussian_filter_length + variables->nrows*variables->ncols;
	complex double* padded_gaussian_filter = malloc(filter_length*sizeof(complex double));
	memset(padded_gaussian_filter, 0, filter_length*sizeof(complex double));
	for(i = 0; i < gaussian_filter_length; i++)
	  padded_gaussian_filter[i] = gaussian_filter[i]/filter_length; // Pre-FFT normalization included.

	fftw_complex* gaussian_filter_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
	fftw_plan gaussian_filter_fft_plan = fftw_plan_dft_1d(filter_length, padded_gaussian_filter, gaussian_filter_fft, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(gaussian_filter_fft_plan);
	fftw_destroy_plan(gaussian_filter_fft_plan);

	complex double* bsig_padded = malloc(filter_length*sizeof(complex double));
	memset(bsig_padded, 0, filter_length*sizeof(complex double));
	memcpy(bsig_padded, bsig, variables->nrows*variables->ncols*sizeof(complex double));

	fftw_complex* bsig_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
	// Pre-FFT normalization.
	for(i = 0; i < filter_length; i++)
	  bsig_padded[i] /= filter_length;

	fftw_plan bsig_fft_plan = fftw_plan_dft_1d(filter_length, bsig_padded, bsig_fft, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(bsig_fft_plan);
	fftw_destroy_plan(bsig_fft_plan);

	fftw_complex* convolved_fft = fftw_malloc(filter_length*sizeof(complex double));

	for(i = 0; i < filter_length; i++)
	  convolved_fft[i] = bsig_fft[i]*gaussian_filter_fft[i];

	fftw_complex* convolved = fftw_malloc(variables->nrows*variables->ncols*sizeof(fftw_complex));
	fftw_plan convolved_ifft_plan = fftw_plan_dft_1d(variables->nrows*variables->ncols, convolved_fft, convolved, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(convolved_ifft_plan);
	fftw_destroy_plan(convolved_ifft_plan);

	// Free allocated in-between memory.
	fftw_free(convolved_fft);
	fftw_free(bsig_fft);
	fftw_free(gaussian_filter_fft);
	free(bsig_padded);
	free(padded_gaussian_filter);
	free(gaussian_filter);
	free(bsig);

	// Pass through integrator -> Output is the cosine argument
	double* phase = malloc(variables->nrows*variables->ncols*sizeof(double));
	for(i = 0; i < variables->ncols*variables->nrows; i++){
	  double sum = 0;
	  int j;
	  for(j = 0; j <= i; j++)
	    sum+= convolved[j];
	  phase[i] = 2*PI*f_c*i+2*PI*f_m*sum;
	}

	// Output is a cosinusoid 
	double rfi_amplitude = 0.1;
	complex double* gsm_signal = malloc(variables->nrows*variables->ncols*sizeof(complex double));
	for(i = 0; i < variables->nrows*variables->ncols; i++){
	  gsm_signal[i] = rfi_amplitude*cos(phase[i]);
	}

	// Add output to radar image
	for(i = 0; i < variables->nrows*variables->ncols; i++){
	  data->radar_image[i] += gsm_signal[i];
	}

	// Free remaining intermediate memory.
	free(phase);
	free(gsm_signal);
}
*/

/* In the processing stages in the simulator, amplitudes may rise dangerously high and cause overflow if not compensated for.
 * This code simply normalizes the input image by the mean.
 */
void normalize_image(complex double* image, unsigned int rows, unsigned int cols){
	double* colmeans = malloc(cols*sizeof(double));
	double rowmean = 0;
	int i,j;
	for(i = 0; i < cols; i++){
	  rowmean = 0;
	  for(j = 0; j < rows; j++){
	    rowmean+= image[i*rows+j];
	  }
	  rowmean/= rows;
	  colmeans[i] = rowmean;
	}
	double mean = 0;
	for(i = 0; i < cols; i++){
	  mean += colmeans[i];
	}

	free(colmeans);

	for(i = 0; i < cols; i++){
	  for(j = 0; j < rows; j++){
	    image[i*rows+j] /= mean;
	  }
	}
}
