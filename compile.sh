#!/bin/bash
gcc waveforms.c algorithms.c sar_simulator.c file_io.c filters.c -lfftw3 -lm -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -march=native -mcmodel=medium -g -o sar_simulator
