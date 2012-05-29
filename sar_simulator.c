/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This simulator can either simulate and process simulated data,
 * or read pre-existing radar data from file, and then process it.
 *
 * The simulator utilizes standard SAR imaging and filtering algorithms, and the
 * processing chain typically looks something like this:
 * create an empty scene -> generate chirp waveform -> generate matched chirp waveform -> 
 * pulse compress single chirp waveform -> insert original (uncompressed) waveform in the middle of the scene ->
 * employ radar scanning algorithm -> employ RFI filtering (if needed) ->
 * employ pulse compression (if needed) -> employ GBP algorithm -> employ apodization (if needed) ->
 * calculate 2D FFT of GBP image -> write *all* data to file
 *
 * This file (sar_simulator.c) contains no algorithms, it just helps the user input data needed for simulation and
 * processing.
 */

#include "sar_simulator.h"

int main(int argc, char** argv){
  radar_variables variables;
  data_arrays data;
 
  memset(&data, 0, sizeof(data_arrays));

  printf("Do you wish to simulate or process radar data? (s/p): ");
  variables.mode = getchar();
  int ret;
  if(variables.mode == 'p'){
    printf("Please enter file name of raw data: ");
    ret = scanf("%s", variables.radar_data_filename);
    if(ret == EOF){
	printf("Invalid input detected, closing.\n");
	return;
    }
    if(!read_radar_file(&data, &variables)){
	printf("Failed to read radar data, closing.\n");
	return;
    }
    process_data(&data, &variables);
  }
  else if(variables.mode == 's'){
    printf("Simulate with real or complex values? (r/c): ");
    ret = scanf("%s", variables.real_or_complex_simulation);
    if(*variables.real_or_complex_simulation != 'r')
      *variables.real_or_complex_simulation = 'c';

    printf("Antenna azimuth beamwidth in radians: ");
    ret = scanf("%f", &variables.beamwidth);

    printf("Chirp start frequency: ");
    ret = scanf("%li", &variables.start_frequency);

    if(variables.start_frequency < 0){
	printf("Negative start frequency entered, closing.\n");
	return;
    }

    printf("Chirp bandwidth: ");
    ret = scanf("%li", &variables.bandwidth);

    if(variables.bandwidth < 0){
	printf("Negative bandwidth entered, closing.\n");
	return;
    }

    printf("Chirp bandwidth-time product: ");
    ret = scanf("%u", &variables.btproduct);
    
    if(variables.btproduct < 1){
	printf("Too small BT-product entered, closing.\n");
	return;
    }

    ret = simulate(&data, &variables);
    if(ret == -1)
      return;
    process_data(&data, &variables);
  }
  else{
    printf("Mode not recognized - exiting.\n");
    return;
  }

  printf("Number of rows in final scene: %i\n", variables.nrows);
  printf("Number of columns in final scene: %i\n", variables.ncols);
  printf("Number of complex points: %i\n",variables.ncols*variables.nrows);
  printf("Uncompressed pulse resolution: %fm\n", variables.signal_distance);

  if(*variables.real_or_complex_simulation == 'c'){
    ret = write_complex_data(&data, &variables);
  }
  else{
    ret = write_real_data(&data, &variables);
  }

}

void free_memory(data_arrays* data){
  if(data->chirp_time_vector)
    free(data->chirp_time_vector);
  if(data->matched_time_vector)
    free(data->matched_time_vector);
  if(data->chirp_fft)
    free(data->chirp_fft);
  if(data->matched_fft);
    free(data->matched_fft);
  if(data->pulse_compressed_waveform)
    free(data->pulse_compressed_waveform);
  if(data->chirp_signal)
    free(data->chirp_signal);
  if(data->matched_chirp)
    free(data->matched_chirp);
  if(data->scene_with_waveform)
    free(data->scene_with_waveform);
  if(data->undistorted_radar_image)
    free(data->undistorted_radar_image);
  if(data->radar_image)
    free(data->radar_image);
  if(data->unfiltered_radar_image)
    free(data->unfiltered_radar_image);
  if(data->pulse_compressed_radar_image)
    free(data->pulse_compressed_radar_image);
  if(data->sar_image)
    free(data->sar_image);
  if(data->sar_fft)
    free(data->sar_fft);
  if(data->sar_img_shifted)
    free(data->sar_img_shifted);
}

int simulate(data_arrays* data, radar_variables* variables){
  chirp_generator(data, variables);
  printf("Chirp generated.\n");

  chirp_matched_generator(data, variables);
  printf("Matched chirp generated.\n");

  data->chirp_fft = malloc(variables->chirp_length*sizeof(complex double));
  fft_waveform(variables->chirp_length, data->chirp_signal, data->chirp_fft);
  printf("Chirp FFT generated.\n");

  data->matched_fft = malloc(variables->chirp_length*sizeof(complex double));
  fft_waveform(variables->chirp_length, data->matched_chirp, data->matched_fft);
  printf("Matched chirp FFT generated.\n");

  pulse_compress_signal(data, variables);
  printf("Pulse-compressed single chirp signal.\n");

  printf("Compressed pulse resolution: %lfm\n", calculate_compressed_pulse_resolution(data, variables));
  
  printf("The target will be placed in the middle of the simulated area.\n");
  printf("Enter area azimuth length (m): ");
  float len = 0;
  int ret = 0;
  ret = scanf("%f", &len);
  variables->ncols = len*variables->chirp_length/variables->signal_distance;
  if(variables->ncols < 2){
    printf("Invalid azimuth length, exiting.\n");
    return -1;
  }
  printf("Enter area range (m): ");
  ret = scanf("%f", &len);
  variables->nrows = len*variables->chirp_length/variables->signal_distance;
  if(variables->nrows < variables->chirp_length){
    printf("Too small range, exiting.\n");
    return -1;
  }

  printf("Please input SAR platform height: ");
  ret = scanf("%d", &variables->altitude);
  if(ret == EOF){
    printf("Invalid input detected, closing.\n");
    return;
  }
  
  printf("Do you want to simulate an object in the scene (y/n)? ");
  char pc = 0;
  do{
    ret = scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);
  if(pc == 'y'){
    printf("Inserting waveform(s) in scene ... ");
    insert_waveform_in_scene(data, variables);
    printf("done.\n");
  }

  printf("Scene range: %fm\n", variables->signal_distance*(variables->nrows/variables->chirp_length));
  printf("Scene azimuth length: %fm\n", variables->signal_distance*(variables->ncols/variables->chirp_length));

  printf("Running radar imager ... ");
  radar_imager(data, variables);
  printf("done.\n");

  printf("Do you want to simulate radio-frequency interference (y/n)? ");
  do{
    ret = scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);
  if(pc == 'y'){
    printf("Adding GSM interference ... ");
    generate_gsm_interference(data, variables);
    printf("done.\n");
  }
}

void process_data(data_arrays* data, radar_variables* variables){
  char pc = 0;
  int ret;

  printf("Do you want to employ RFI suppression (y/n)? ");
  do{
    ret = scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);
  if(pc == 'y'){
    printf("Suppressing RFI ... ");
    rfi_suppression(data, variables);
    printf("done.\n");
  }
  else if(pc == 'n'){
    if(!data->unfiltered_radar_image){
	data->unfiltered_radar_image = malloc(variables->ncols*variables->nrows*sizeof(complex double));
	memcpy(data->unfiltered_radar_image, data->radar_image, variables->ncols*variables->nrows*sizeof(complex double));
    }
  }
  
  printf("Do you want to enable pulse compression (y/n)? ");
  do{
    ret = scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);
  if(pc == 'y'){
    printf("Pulse-compressing image ... ");
    pulse_compress_image(data, variables);
    printf("done.\n");
  }
  else if(pc == 'n'){
    data->pulse_compressed_radar_image = malloc(variables->nrows*variables->ncols*sizeof(complex double));
    memcpy(data->pulse_compressed_radar_image, data->radar_image, variables->nrows*variables->ncols*sizeof(complex double));
  }
  
  printf("Running GBP ... ");
  gbp(data, variables);
  printf("done.\n");

  printf("Generating 2D FFT of GBP image ... ");
  gbp_fft(data, variables);
  printf("done.\n");

  if(variables->mode == 's')
    filter_dc(data, variables);

  /*
  printf("Do you want to apodize SAR image (y/n)? ");
  do{
    ret = scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);

  if(pc == 'y'){
    printf("Apodizing SAR image ... ");
    apodize_sar_fft(data, variables);
    printf("done.\n");
  }
  */
}
