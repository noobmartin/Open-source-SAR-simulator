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

void write_data(matrix* data_matrix){
  FILE* fp;
  FILE* dimensions;
  
  char filename[255];
  dimensions = fopen("output/dimensions.dat", "a");

  memset(filename, 0, 255);
  strcat(filename, data_matrix->name);
  strcat(filename, ".dat");
  
  fp = fopen(filename, "w");

  for(int i = 0; i < data_matrix->cols; i++){
    for(int j = 0; j < data_matrix->rows; j++){
      fprintf(fp, FILE_OUTPUT_PRECISION, creal(data_matrix->data[i*data_matrix->rows+j]));
      fprintf(fp, FILE_OUTPUT_PRECISION, cimag(data_matrix->data[i*data_matrix->rows+j]));
    }
    fprintf(fp, "\n");
  }
  
  fclose(fp);

  fprintf(dimensions, "%s\n%i\n%i\n", data_matrix->name, data_matrix->rows, data_matrix->cols);
  
  fclose(dimensions);

  return;
}
