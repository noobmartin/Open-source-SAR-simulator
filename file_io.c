/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This code is ugly and really ought to be re-written.
 */

#include "sar_simulator.h"

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
 
  int i,j;
  FILE* fp;
  FILE* dimensions;
  matrix* data_ptr = data_matrix;
  char filename[255];
  dimensions = fopen("dimensions.dat", "w");
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
      for(i = 0; i < data_ptr->cols; i++){
	for(j = 0; j < data_ptr->rows; j++){
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

  //fseek(fp, 0L, SEEK_END);
  //unsigned int file_size = ftell(fp);
  //fseek(fp, 0L, SEEK_SET);
  
  radar_metadata meta;
  
  int ret;
  printf("Radar rows: ");
  ret = scanf("%u", &meta.rows);

  printf("Radar cols: ");
  ret = scanf("%u", &meta.cols);

  meta.real_or_complex = 'r';

  matrix* meta_radar = get_last_node(data);
  strcpy(meta_radar->name, "radar_image");
  meta_radar->rows = meta.rows;
  meta_radar->cols = meta.cols;
  meta_radar->data = malloc(meta.rows*meta.cols*sizeof(complex double));
  memset(meta_radar->data, 0, meta.rows*meta.cols*sizeof(complex double));
  complex double* radar_image = meta_radar->data;

  //FILE* mp = fopen("radar_metadata", "r");
  //ret = fread(&meta, sizeof(radar_metadata), 1, mp);

  //fclose(mp);


  int i,j;
  unsigned int real;
  complex double imag;
  if(meta.real_or_complex == 'r'){
      for(j = 0; j < meta.cols; j++){
    for(i = 0; i < meta.rows; i++){
        //ret = fread(data->radar_image + i*sizeof(complex double), sizeof(double), 1, fp);
	ret = fscanf(fp, "%u", &real);
	imag = real + _Complex_I*0;
	radar_image[j*meta.rows+i] = imag;
      }
    }
  }
  else if(meta.real_or_complex == 'c'){
      for(j = 0; j < meta.rows; j++){
    for(i = 0; i < meta.rows*meta.cols; i++){
      //ret = fread(data->radar_image + i*sizeof(complex double), sizeof(complex double), 1, fp);
        ret = fscanf(fp, "%lf", (double*)(  &radar_image[i*meta.rows+j]  ));
      }
    }
  }
  else{
    printf("Invalid data mode, should be real or complex, read %c.\n", meta.real_or_complex);
    fclose(fp);
    return -1;
  }

  fclose(fp);
  return 1;
}
