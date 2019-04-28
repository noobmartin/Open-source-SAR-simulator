/* Author: Alexander Rajula
 * Contact: alexander@rajula.org
 *
 * This header file contains basic structures for managing all data structures in the simulator.
 */

#ifndef sar_simulator_h
#define sar_simulator_h

#include "common_types.h"

void simulate(matrix* data_matrix, radar_variables* variables);
void process_data(matrix* data_matrix, radar_variables* variables);
void free_memory(matrix* data_matrix);
void build_metadata(matrix* data);

#endif
