/* This file contains filters extracted from the CinSnow
 * radar application.
 * Credits are due to Dan Axelsson @ Cinside and M. Larsson
 */

#include "cinsnowfilters.h"
#include "sar_simulator.h"

double *gpAvg;
double *gpAvgBG;
double* IIR_filtOut;
double iirFilterConst;
double lpFilter;
unsigned int bgWindowSize;
unsigned int rollingAverageSize;

void cinsnowfilters(matrix* data, radar_variables* variables){
 int ret;
 //printf("Enter IIR filter constant: ");
 //ret = scanf("%lf", &iirFilterConst);
 printf("Enter LP filter constant: ");
 ret = scanf("%lf", &lpFilter);
 printf("Enter rolling average constant: ");
 ret = scanf("%u", &rollingAverageSize);
 printf("Enter background subtraction constant: ");
 ret = scanf("%u", &bgWindowSize);

  matrix* meta_radar = get_matrix(data, "radar_image");
  complex double* radar_image = meta_radar->data;

 int i;
 double* column;
 for(i = 0; i < meta_radar->cols; i++){
    column = (double*)&radar_image[i*meta_radar->rows];
    
    //IIR_filter(column, variables->nrows);
    subtract_dc_offset(column, meta_radar->rows*2);
    low_pass_filter(column, meta_radar->rows*2, lpFilter);
    //DSubref(column, gTraceRef, variables->nrows);
    rolling_average(column, meta_radar->rows*2, rollingAverageSize);
    update_background_map(column, meta_radar->rows*2, bgWindowSize);
    DSubref(column, gpAvgBG, meta_radar->rows*2);
 }
 
 int j;
 for(i = 0; i < 80; i++){
 for(j = 0; j < 35; j++){
   column = (double*)&radar_image[j*meta_radar->rows];
   update_background_map(column, meta_radar->rows*2, bgWindowSize);
   DSubref(column, gpAvgBG, meta_radar->rows*2);
 }
 }
 
 normalize_image(radar_image, meta_radar->rows, meta_radar->cols);
}

void DSubref(double* data, double* ref, unsigned int number){
  int i;
  if(ref == NULL)
    return;
  for(i = 0; i < number; i++)
    *data++ -= *ref++;
}

void IIR_filter(double* traceData, unsigned int samples){
  int i;
  static char bInit = 0;

  // Initialize filter.
  if(bInit == 0){
    IIR_filtOut = malloc(samples*sizeof(double));
    for(i = 0; i < samples; i++){
      IIR_filtOut[i] = traceData[i];
    }
    bInit = 1;
  }
  // Update filter.
  else{
    for(i = 0; i < samples; i++){
      IIR_filtOut[i] = IIR_filtOut[i]*(1-iirFilterConst)+traceData[i]*iirFilterConst;
      traceData[i] = (abs(traceData[i] - IIR_filtOut[i]));
    }
  }
}

void subtract_dc_offset(double* data, unsigned int size){
  int i, DCoffset;
  DCoffset = 0;
  for(i = 0; i < size; i++){
    DCoffset+= data[i];
  }
  DCoffset/=size;
  for(i = 0; i < size; i++)
    data[i]-= DCoffset;
}

void low_pass_filter(double* data, unsigned int size, int filter){
  int i,d;
  d = data[0];
  for(i = 0; i < size; i++){
    d = (int)((double)(filter*d + (long)(16-filter)*data[i] + 8)/16);
    data[i] = d;
  }
}

void rolling_average(double* data, unsigned int size, int mavg){
  int i;
  if(gpAvg == NULL){
    gpAvg = malloc(size*sizeof(double));
    memset(gpAvg, 0, size*sizeof(double));
  }

  for(i = 0; i < size; i++){
    data[i] = gpAvg[i] = (int)((mavg*(int)gpAvg[i] + (16-mavg)*(int)data[i] + 8)/16);
  }
}

void update_background_map(double* data, unsigned int size, int mavg){
  int i;
  if(gpAvgBG == NULL){
    gpAvgBG = malloc(size*sizeof(double));
    memset(gpAvgBG, 0, size*sizeof(double));
    for(i = 0; i < size; i++)
      gpAvgBG[i] = data[i];
  }

  for(i = 0; i < size; i++){
    gpAvgBG[i] = (int)(1.0/mavg*data[i] + (mavg-1.0)/mavg*gpAvgBG[i]);
  }
}
