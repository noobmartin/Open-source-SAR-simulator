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
 */
void gbp(matrix* data, radar_variables* variables){
  matrix* meta_radar_image;

  if(variables->mode == Process)
    meta_radar_image = get_matrix(data, "radar_image");
  else
    meta_radar_image = get_matrix(data, "pc_image");

  matrix* meta_sar_image = get_last_node(data);
  strcpy(meta_sar_image->name, "sar_image");

  complex double* radar_image = meta_radar_image->data;

	meta_sar_image->rows = meta_radar_image->rows;
	meta_sar_image->cols = meta_radar_image->cols;
	meta_sar_image->data = malloc(meta_sar_image->rows*meta_sar_image->cols*sizeof(double complex));

  unsigned int range_index;
	unsigned int cols = meta_sar_image->cols;
	unsigned int rows = meta_sar_image->rows;

  printf("Running GBP.\n");

  for(unsigned int j = 0; j < cols; j++){
    printf("Pass %i of %i.\n", j, cols);

    for(unsigned int k = 0; k < rows; k++){
      for(unsigned int l = 0; l < cols; l++){
				range_index = sqrt((l-j)*(l-j)+k*k);
				if(range_index < rows){
					meta_sar_image->data[j*rows+k] += radar_image[l*rows+range_index];
				}
      }
    }
	}

	printf("GBP done.\n");
}

/* To evaluate the quality of a simulated SAR image, one often looks at the spectrum of the GBP image.
 * For a point object in the GBP image, the 2D FFT should have a "fan" shape. */
void gbp_fft(matrix* data, radar_variables* variables){
  matrix* meta_sar_image = get_matrix(data, "sar_image");

  unsigned int rows = meta_sar_image->rows;
  unsigned int cols = meta_sar_image->cols;

  matrix* meta_sar_fft = get_last_node(data);
  strcpy(meta_sar_fft->name, "sar_fft");
  
  meta_sar_fft->rows = rows;
  meta_sar_fft->cols = cols;
  meta_sar_fft->data = malloc(rows*cols*sizeof(double complex));
  
  complex double* sar_img_shifted = malloc(rows*cols*sizeof(double complex));
  complex double* sar_image = meta_sar_image->data;

  // This shifts the frequencies of the 2D FFT so that zero frequency is in the middle of the image.
  // Also pre-normalizes the data since the 2D FFT scales by a factor of rows*cols.
  for(unsigned int i = 0; i < cols; i++){
    for(unsigned int j = 0; j < rows; j++){
			sar_img_shifted[rows*i+j] = sar_image[rows*i+j]*pow(-1,i+j)/(cols*rows);
    }
  }
  
  fftw_plan fft = fftw_plan_dft_2d(cols, rows, sar_img_shifted, meta_sar_fft->data, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);

  free(sar_img_shifted);
}

/* This radar imaging algorithm is designed by me, and simulates scanning over a scene.
 * It takes antenna azimuth beamwidth and antenna elevation into account, and calculates how many azimuth bins the antenna sees from each azimuth position.
 * All the data within the antenna's field of view is compressed into a single column, for each column, taking distance into account.
 * This is what generates the familiar hyperbola in the raw radar image.
 */
void radar_imager(matrix* data, radar_variables* variables){
  printf("Which algorithm should I use to simulate the scan?\n");
  printf("1) Radar imaging algorithm\n");
  printf("2) Coherent standard algorithm, time delay method\n");
  printf("3) Coherent standard algorithm, phase shift method\n");
  
  unsigned int algorithm = 0;
  scanf("%u", &algorithm);

  matrix* meta_scene = get_matrix(data, "scene");

  double complex* scene = meta_scene->data;

  unsigned int rows = meta_scene->rows;
  unsigned int cols = meta_scene->cols;
 
  matrix* meta_radar_image = get_last_node(data);
  strcpy(meta_radar_image->name, "radar_image");
  meta_radar_image->rows = rows;
  meta_radar_image->cols = cols;
  meta_radar_image->data = malloc(rows*cols*sizeof(complex double));
  double complex* radar_image = meta_radar_image->data;

  matrix* meta_chirp = get_matrix(data, "chirp");

  complex double* chirp = meta_chirp->data;
  
  if(algorithm == 1){
  	printf("Using Radar imaging algorithm.\n");

	  printf("Enter antenna azimuth beamwidth in radians: ");
	  scanf("%f", &variables->beamwidth);

	  printf("Enter SAR platform height: ");
	  scanf("%d", &variables->altitude);

	  unsigned int i,j,k;
	  int beam_value;
	  double azimuth_coverage = round(variables->altitude*tan(0.5*variables->beamwidth));
	  unsigned int beamcrossrange = variables->chirp_length*azimuth_coverage/variables->signal_distance;
	  
	  printf("Antenna beam covers %f meters, %i columns.\n", azimuth_coverage, beamcrossrange);
	  for(i = 0; i < cols; i++){
	    printf("Executing RIA, pass %i of %i.\n", i+1, cols);
	    
	    for(j = 0; j < 2*beamcrossrange; j++){
	      beam_value = j-beamcrossrange;
	      for(k = 0; k < rows; k++){
		      if(i+beam_value < cols){
		        unsigned int dist = sqrt(pow(beam_value,2)+pow(k, 2));
		        if(dist < rows){
		          if(i+beam_value >= 0){
		            radar_image[i*rows+dist] += scene[(i+beam_value)*rows+k];
		          }
		        }
		      }
	      }
	    }
	  }

	  printf("Radar imaging algorithm done.\n");
  }
  else if(algorithm == 2){
  	printf("Using coherent standard algorithm, method one.\n");

	  for(unsigned int i = 0; i < cols; i++){
	    unsigned int dist = 0;
	    
	    if(i < cols/2){
	      dist = sqrt(pow(cols/2-i,2)+pow(rows/2-variables->chirp_length/2,2));
	    }
	    else{
	      dist = sqrt(pow(i-cols/2,2)+pow(rows/2-variables->chirp_length/2,2));
	    }

	    if(dist+variables->chirp_length < rows){
	      memcpy(&radar_image[i*rows+dist], chirp, meta_chirp->rows*sizeof(complex double));
	      
	      for(unsigned int j = 0; j < meta_chirp->rows; j++){
	        radar_image[i*rows+dist+j] /= pow(dist,4);
	      }
	    }
	  }
	  printf("Standard algorithm done.\n");
	}
	else{
 	  printf("Using coherent standard algorithm, method two.\n");

	  memset(radar_image, 0, rows*cols*sizeof(complex double));

	  // Insert sent waveform.
	  for(unsigned int i = 0; i < cols; i++){
	    memcpy(&radar_image[i*rows], chirp, meta_chirp->rows*sizeof(complex double));
	  }

	  double sample_distance = variables->signal_distance/variables->chirp_length;
	
	  double* freq      = malloc(rows*sizeof(double));	
	  double freq_delta = 5*variables->bandwidth/rows;
	  double freq_inc   = 0;

	  for(unsigned int i = 0; i < rows; i++){
	    freq[i] = freq_inc;
	    freq_inc += freq_delta;
	  }

	  // Create the hyperbola.
	  for(unsigned int i = 0; i < cols; i++){
	    unsigned int dist = 0;
	    if(i < cols/2)
	      dist = sqrt(pow(cols/2-i,2)+pow(rows/2-meta_chirp->rows/2,2));
	    else
	      dist = sqrt(pow(i-cols/2,2)+pow(rows/2-meta_chirp->rows/2,2));

	    fftw_plan fft = fftw_plan_dft_1d(rows, &radar_image[i*rows], &radar_image[i*rows], FFTW_FORWARD, FFTW_ESTIMATE);
	    fftw_execute(fft);
	    fftw_destroy_plan(fft);

	    double in_transit_time = dist*sample_distance/C;
	    for(unsigned int j = 0; j < rows; j++){
	      // The division of "rows" here is to normalize the FFT->IFFT operation, since the FFTW library doesn't re-normalize by default.
	      radar_image[i*rows+j] *= cexp(-_Complex_I*2*PI*freq[j]*in_transit_time)/(rows*pow(dist,4));
	    }
	    
	    fftw_plan ifft = fftw_plan_dft_1d(rows, &radar_image[i*rows], &radar_image[i*rows], FFTW_BACKWARD, FFTW_ESTIMATE);
	    fftw_execute(ifft);
	    fftw_destroy_plan(ifft);

	  }

	  free(freq);

	  printf("Standard algorithm done.\n");
   }

}

// This just copies the waveform into the center of the scene, to simulate an object at the center.
void insert_waveform_in_scene(matrix* data, radar_variables* variables){
  printf("The target will be placed in the middle of the simulated area.\n");
  
  matrix* meta_chirp = get_matrix(data, "chirp");

  matrix* meta_scene = get_matrix(data, "scene");
  if(meta_scene == NULL){
    meta_scene = get_last_node(data);
    strcpy(meta_scene->name, "scene");
  }

  printf("Enter area azimuth length (m): ");
  double len = 0;
  scanf("%lf", &len);
  meta_scene->cols = len*meta_chirp->rows/variables->signal_distance;
  
  printf("Enter area range (m): ");
  scanf("%lf", &len);
  meta_scene->rows = len*meta_chirp->rows/variables->signal_distance;
  
  if(meta_scene->cols < 2){
    printf("Invalid azimuth length, exiting.\n");
    return;
  }
  
  if(meta_scene->rows < meta_chirp->rows){
    printf("Too small range, exiting.\n");
    return;
  }

  unsigned int rows = meta_scene->rows;
  unsigned int cols = meta_scene->cols;
  
  meta_scene->data = malloc(rows*cols*sizeof(complex double));
  complex double* scene = meta_scene->data;
 
  printf("Scene range: %fm\n",          variables->signal_distance*(rows/meta_chirp->rows));
  printf("Scene azimuth length: %fm\n", variables->signal_distance*(cols/meta_chirp->rows));
 
  memcpy(&scene[(unsigned int)(cols/2)*rows + (unsigned int)(rows/2) - meta_chirp->rows/2], 
         meta_chirp->data, 
         meta_chirp->rows*sizeof(complex double));
}

/* Pulse compression of a single waveform (it's nice to look at the plots to see that pulse compression works).
 * This is just convolution computed in the frequency domain with the FFT.
 */
void pulse_compress_signal(matrix* data, radar_variables* variables){
  matrix* meta_pc_waveform = get_last_node(data);
  matrix* meta_chirp       = get_matrix(data, "chirp");
  matrix* meta_match       = get_matrix(data, "match");
  
  strcpy(meta_pc_waveform->name, "pulse_compressed_waveform");
  meta_pc_waveform->rows = meta_chirp->rows;
  meta_pc_waveform->cols = meta_chirp->cols;
  meta_pc_waveform->data = malloc(meta_chirp->rows*sizeof(complex double));
  
  unsigned int kernel_length = variables->chirp_length;
  unsigned int filter_length = 2*kernel_length;

  fftw_complex* padded_signal = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  
  fftw_complex* sigfft   = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* matchfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product  = fftw_malloc(filter_length*sizeof(fftw_complex));
  
  fftw_plan sigp   = fftw_plan_dft_1d(filter_length, padded_signal, sigfft,                 FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan matchp = fftw_plan_dft_1d(filter_length, padded_kernel, matchfft,               FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan iff    = fftw_plan_dft_1d(kernel_length, product,       meta_pc_waveform->data, FFTW_FORWARD, FFTW_MEASURE);
  
  memset(padded_signal, 0, filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_signal, meta_chirp->data, kernel_length*sizeof(fftw_complex));
  memcpy(padded_kernel, meta_match->data, kernel_length*sizeof(fftw_complex));

  for(unsigned int i = 0; i < filter_length; i++){
    padded_signal[i] /= filter_length;
    padded_kernel[i] /= filter_length;
  }

  fftw_execute(sigp);
  fftw_execute(matchp);

  for(unsigned int i = 0; i < filter_length; i++){
    product[i] = sigfft[i]*pow(-1,i)*conj(sigfft[i]);
  }

  fftw_execute(iff);

  fftw_destroy_plan(sigp);
  fftw_destroy_plan(matchp);
  fftw_destroy_plan(iff);

  fftw_free(sigfft);
  fftw_free(matchfft);
  fftw_free(product);

}

/* Convolution of the entire radar image with the matched waveform, again implemented with the FFT. */
void pulse_compress_image(matrix* data, radar_variables* variables){
  matrix* meta_radar_image = get_matrix(data, "radar_image");
  complex double* radar_image = meta_radar_image->data;
  
  matrix* meta_pc_image = get_last_node(data);
  strcpy(meta_pc_image->name, "pc_image");

  matrix* meta_match = get_matrix(data, "match");
  complex double* match = meta_match->data;

  unsigned int rows = meta_radar_image->rows;
  unsigned int cols = meta_radar_image->cols;
  
  meta_pc_image->rows = rows;
  meta_pc_image->cols = cols;
  
  meta_pc_image->data = malloc(rows*cols*sizeof(double complex));
  complex double* pc_image = meta_pc_image->data;

  // Make sure input has valid values.
  unsigned int kernel_length = meta_match->rows;
  unsigned int z;

  normalize_image(radar_image, rows, cols);

  unsigned int filter_length = rows + kernel_length;

  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* kernel_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memset(kernel_fft, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_kernel, match, kernel_length*sizeof(fftw_complex));

  // Normalize before FFT.
  for(z = 0; z < kernel_length; z++){
    padded_kernel[z] /= filter_length;
  }

  // Compute fft of filter kernel.
  fftw_plan kernelfft = fftw_plan_dft_1d(filter_length, padded_kernel, kernel_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(kernelfft);
  fftw_complex* output_column;
  fftw_complex* padded_column = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_column_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  unsigned int i,j;
  for(i = 0; i < cols; i++){
    double complex* column = &radar_image[i*rows];

    output_column = &pc_image[i*rows];
    memset(padded_column, 0, filter_length*sizeof(fftw_complex));
    memcpy(padded_column, column, rows*sizeof(fftw_complex));

    // Normalize before FFT.
    for(z = 0; z < filter_length; z++){
      padded_column[z] /= filter_length;
    }

    memset(padded_column_fft, 0, filter_length*sizeof(fftw_complex));
    fftw_plan colfft = fftw_plan_dft_1d(filter_length, padded_column, padded_column_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(colfft);

    for(j = 0; j < filter_length; j++){
			product[j] = padded_column_fft[j]*kernel_fft[j];
    }

    fftw_plan colifft = fftw_plan_dft_1d(rows, product, output_column, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(colifft);

    fftw_destroy_plan(colfft);
    fftw_destroy_plan(colifft);
  }

  fftw_free(product);
  fftw_free(padded_column_fft);
  fftw_free(padded_column);
  fftw_free(padded_kernel);

  fftw_destroy_plan(kernelfft);

}

void pulse_compress_scene(matrix* data, radar_variables* variables){
  matrix* meta_match = get_matrix(data, "match");
  matrix* meta_scene = get_matrix(data, "scene");

  unsigned int rows = meta_scene->rows;
  unsigned int cols = meta_scene->cols;

  matrix* meta_pc_scene = get_last_node(data);
  strcpy(meta_pc_scene->name, "pc_scene");
  meta_pc_scene->rows = rows;
  meta_pc_scene->cols = cols;
  meta_pc_scene->data = malloc(cols*rows*sizeof(double complex));
  
  unsigned int kernel_length = meta_match->rows;
  unsigned int filter_length = rows + kernel_length;

  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* kernel_fft    = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product       = fftw_malloc(filter_length*sizeof(fftw_complex));
  
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memset(kernel_fft, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_kernel, meta_match->data, kernel_length*sizeof(fftw_complex));

  // Compute fft of filter kernel.
  fftw_plan kernelfft = fftw_plan_dft_1d(filter_length, padded_kernel, kernel_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(kernelfft);
  fftw_complex* output_column;
  fftw_complex* padded_column = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_column_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  
  for(unsigned int i = 0; i < cols; i++){
    double complex* column = &meta_scene->data[i*rows];

    output_column = &meta_pc_scene->data[i*rows];
    memset(padded_column    , 0     , filter_length*sizeof(fftw_complex));
    memcpy(padded_column    , column, rows*sizeof(fftw_complex));
    memset(padded_column_fft, 0     , filter_length*sizeof(fftw_complex));
    
    fftw_plan colfft = fftw_plan_dft_1d(filter_length, padded_column, padded_column_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(colfft);

    for(unsigned int j = 0; j < filter_length; j++){
			product[j] = padded_column_fft[j]*kernel_fft[j];
    }

    fftw_plan colifft = fftw_plan_dft_1d(rows, product, output_column, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(colifft);

    fftw_destroy_plan(colfft);
    fftw_destroy_plan(colifft);
  }

  fftw_free(product);
  fftw_free(padded_column_fft);
  fftw_free(padded_column);
  fftw_free(padded_kernel);

  fftw_destroy_plan(kernelfft);

  memcpy(meta_scene->data, meta_pc_scene->data, rows*cols*sizeof(double complex));
}
