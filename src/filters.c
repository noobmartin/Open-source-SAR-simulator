/*
 * Author: Alexander Rajula
 * Contact: alexander@rajula.org
 */

#include "filters.h"
#include <stdlib.h>

/* In the processing stages in the simulator, amplitudes may rise dangerously high and cause overflow if not compensated for.
 * This code simply normalizes the input image by the mean.  */
void normalize_image(complex double* image, unsigned int rows, unsigned int cols){
	double* colmeans = malloc(cols*sizeof(complex double));
	double  rowmean  = 0;
	double  mean     = 0;

	for(int i = 0; i < cols; i++){
	  rowmean = 0;

	  for(int j = 0; j < rows; j++){
	    rowmean+= image[i*rows+j];
	  }

	  rowmean /= rows;
	  colmeans[i] = rowmean;
	}

	for(int i = 0; i < cols; i++){
	  mean += colmeans[i];
	}
	
	mean /= cols;

	free(colmeans);

	for(int i = 0; i < cols; i++){
	  for(int j = 0; j < rows; j++){
	    image[i*rows+j] /= mean;
	  }
	}

}
