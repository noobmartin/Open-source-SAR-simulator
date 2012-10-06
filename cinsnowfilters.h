#ifndef _cinsnowfilters_h
#define _cinsnowfilters_h

#include "sar_simulator.h"

void cinsnowfilters(matrix* data, radar_variables* variables);
void DSubref(double* data, double* ref, unsigned int number);
void IIR_filter(double* traceData, unsigned int samples);
void subtract_dc_offset(double* data, unsigned int size);
void low_pass_filter(double* data, unsigned int size, int filter);
void rolling_average(double* data, unsigned int size, int mavg);
void update_background_map(double* data, unsigned int size, int mavg);

#endif
