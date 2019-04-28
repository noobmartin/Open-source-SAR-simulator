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

int main(int argc, char** argv){
  radar_variables variables;
  matrix* data = malloc(sizeof(matrix)); 
  memset(data, 0, sizeof(matrix));
  strcpy(data->name, "metadata");

  printf("Do you wish to simulate or process radar data? (s/p): ");
  
  if(getchar() == 'p'){
    variables.mode = Process;
    
    printf("Please enter file name of raw data: ");

    if(scanf("%s", variables.radar_data_filename) == EOF){
			printf("Invalid input detected, closing.\n");
			return 0;
    }

    if(!read_radar_file(data, &variables)){
			printf("Failed to read radar data, closing.\n");
			return 0;
    }

  }
  else{
    variables.mode = Simulate;
    
    simulate(data, &variables);
  }

  process_data(data, &variables);
  
  build_metadata(data);

  write_data(data, &variables);
 
  //free_memory(data);

  return 0;
}

void free_memory(matrix* data_matrix){
  matrix* ptr = data_matrix;
  matrix* next = ptr->next;
  while(ptr != NULL){
    if(ptr->data != NULL){
      free(ptr->data);
    }

		next = ptr->next;
		free(ptr);
		ptr = next;
  }
}

void simulate(matrix* data, radar_variables* variables){
  chirp_generator(data, variables);
  chirp_matched_generator(data, variables);
 
  matrix* meta_chirp = get_matrix(data, "chirp");
  matrix* meta_match = get_matrix(data, "match");
  
  matrix* meta_chirp_fft = get_last_node(data);
  meta_chirp_fft->data   = malloc(meta_chirp->rows*sizeof(complex double));
  meta_chirp_fft->rows   = meta_chirp->rows;
  meta_chirp_fft->cols   = meta_chirp->cols;
  strcpy(meta_chirp_fft->name, "chirp_fft");

  matrix* meta_match_fft = get_last_node(data);
  meta_match_fft->data   = malloc(meta_match->rows*sizeof(complex double));
  meta_match_fft->rows   = meta_chirp->rows;
  meta_match_fft->cols   = meta_chirp->cols;
  strcpy(meta_match_fft->name, "match_fft");

  fft_waveform(meta_chirp->rows, meta_chirp->data, meta_chirp_fft->data);
  fft_waveform(meta_match->rows, meta_match->data, meta_chirp_fft->data);

  pulse_compress_signal(data, variables);
  
  printf("Compressed pulse resolution: %lfm\n", calculate_compressed_pulse_resolution(data, variables));

  insert_waveform_in_scene(data,variables);

  radar_imager(data, variables);

  return;
}

void process_data(matrix* data, radar_variables* variables){
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
    pulse_compress_image(data, variables);
    printf("done.\n");
  }
  
  gbp(data, variables);

  printf("Generating 2D FFT of GBP image ... ");
  gbp_fft(data, variables);
  printf("done.\n");
}

void build_metadata(matrix* data){
  matrix* ptr = data->next;
  unsigned int elements = 0;
  while(ptr != NULL){
    elements++;
    ptr = ptr->next;
  }
  data->data = malloc(elements*sizeof(complex double));
}
