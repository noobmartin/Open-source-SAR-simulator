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

#include "common_functions.h"
#include "sar_simulator.h"
#include "file_io.h"
#include "waveforms.h"
#include "algorithms.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

matrix time_vector;
matrix chirp;
matrix match;
matrix pc_waveform;
matrix scene;
matrix radar_image;
matrix pc_image;
matrix sar_image;
matrix sar_image_fft;
matrix chirp_fft;
matrix match_fft;
radar_variables variables;
  
int main(int argc, char** argv){
  strcpy(time_vector.name, "time_vector");
  strcpy(chirp.name, "chirp");
  strcpy(chirp_fft.name, "chirp_fft");
  strcpy(match.name, "match");
  strcpy(match_fft.name, "match_fft");
  strcpy(radar_image.name, "radar_image");
  strcpy(pc_image.name, "pc_image");
  strcpy(sar_image.name, "sar_image");
  strcpy(sar_image_fft.name, "sar_fft");
  strcpy(scene.name, "scene");
  strcpy(pc_waveform.name, "pulse_compressed_waveform");
  
  simulate();

  process_data();
  
  write_data(&time_vector);
  write_data(&chirp);
  write_data(&match);
  write_data(&pc_waveform);
  write_data(&scene);
  write_data(&radar_image);
  write_data(&pc_image);
  write_data(&sar_image);
  write_data(&sar_image_fft);
  write_data(&chirp_fft);
  write_data(&match_fft);
 
  return 0;
}

void simulate(void){
  chirp_generator(&time_vector, &chirp, &variables);
  chirp_matched_generator(&chirp, &match);
  
  chirp_fft.data   = malloc(chirp.rows*sizeof(complex double));
  chirp_fft.rows   = chirp.rows;
  chirp_fft.cols   = chirp.cols;


  match_fft.data   = malloc(match.rows*sizeof(complex double));
  match_fft.rows   = chirp.rows;
  match_fft.cols   = chirp.cols;

  fft_waveform(chirp.rows, chirp.data, chirp_fft.data);
  fft_waveform(match.rows, match.data, chirp_fft.data);

  pulse_compress_signal(&chirp, &match, &pc_waveform, &variables);
  
  printf("Compressed pulse resolution: %lfm\n", calculate_compressed_pulse_resolution(&pc_waveform, &variables));

  insert_waveform_in_scene(&chirp, &scene, &variables);

  radar_imager(&scene, &radar_image, &chirp, &variables);

  return;
}

void process_data(void){
  char pc = 0;

  printf("Do you want to enable pulse compression (y/n)? ");
  do{
    scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);

  if(pc == 'y'){
    printf("Pulse-compressing image ... ");
    pulse_compress_image(&radar_image, &pc_image, &match, &variables);
    printf("done.\n");
  }
  
  gbp(&pc_image, &sar_image, &variables);

  printf("Generating 2D FFT of GBP image ... ");
  gbp_fft(&sar_image, &sar_image_fft, &variables);
  printf("done.\n");
}
