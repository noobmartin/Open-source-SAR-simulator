/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This code is ugly and really ought to be re-written.
 */

#include "file_io.h"
#include "common_functions.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int write_data(matrix* data_matrix, radar_variables* variables){
  char fmode = 0;
  int ret = 0;

  do{
    fmode = getchar();
  }while(fmode != '\n');

  printf("Would you like to write data in human-readable or binary format (h/b): ");
  do{
    ret = scanf("%c", &fmode);
  }while((fmode != 'h') && (fmode != 'b'));
 
  FILE* fp;
  FILE* dimensions;
  matrix* data_ptr = data_matrix;
  
  char filename[255];
  dimensions = fopen("output/dimensions.dat", "w");
  while(data_ptr != NULL){
    memset(filename, 0, 255);
    strcat(filename, data_ptr->name);
    strcat(filename, ".dat");
    fp = fopen(filename, "w");
    if(fp == NULL)
      break;
    if(fmode == 'b'){
      ret = fwrite(data_ptr->data, 1, data_ptr->rows*data_ptr->cols*sizeof(complex double), fp);
    }
    else{
      for(int i = 0; i < data_ptr->cols; i++){
        for(int j = 0; j < data_ptr->rows; j++){
          fprintf(fp, FILE_OUTPUT_PRECISION, creal(data_ptr->data[i*data_ptr->rows+j]));
          fprintf(fp, FILE_OUTPUT_PRECISION, cimag(data_ptr->data[i*data_ptr->rows+j]));
        }
        fprintf(fp, "\n");
      }
    }
    
    fclose(fp);

    fprintf(dimensions, "%s\n%i\n%i\n", data_ptr->name, data_ptr->rows, data_ptr->cols);

    data_ptr = data_ptr->next;
  }

  return 0;
}

int read_radar_file(matrix* data, radar_variables* variables){
  FILE* fp = fopen(variables->radar_data_filename, "r");
  if(fp == NULL)
    return -1;

  radar_metadata meta;
  
  int ret;
  printf("Radar rows: ");
  ret = scanf("%u", &meta.rows);

  printf("Radar cols: ");
  ret = scanf("%u", &meta.cols);

  meta.is_complex = true;

  matrix* meta_radar = get_last_node(data);
  strcpy(meta_radar->name, "radar_image");
  meta_radar->rows = meta.rows;
  meta_radar->cols = meta.cols;
  meta_radar->data = malloc(meta.rows*meta.cols*sizeof(complex double));
  memset(meta_radar->data, 0, meta.rows*meta.cols*sizeof(complex double));
  complex double* radar_image = meta_radar->data;

  unsigned int real;
  complex double imag;
  if(!meta.is_complex){
    for(int j = 0; j < meta.cols; j++){
      for(int i = 0; i < meta.rows; i++){
			  ret = fscanf(fp, "%u", &real);
			  imag = real + _Complex_I*0;
			  radar_image[j*meta.rows+i] = imag;
      }
    }
  }
  else{
    for(int j = 0; j < meta.rows; j++){
      for(int i = 0; i < meta.rows*meta.cols; i++){
        ret = fscanf(fp, "%lf", (double*)(  &radar_image[i*meta.rows+j]  ));
      }
    }
  }

  fclose(fp);
  return 1;
}
