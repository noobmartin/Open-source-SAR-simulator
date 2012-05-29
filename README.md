# Contact
You may contact me (Alexander Rajula) at alexander@rajula.org or superrajula@gmail.com

or by phone (+46)705299713

If you have any questions about the code, or have found it useful,
don't hesitate to e-mail me.

# Copyright
The code and all files within this folder is supplied "as-is".
You are allowed to make copies of it and modify it to your wishes as long as you leave a note in 
your code (or other documentation) about its original author.

# Introduction
The code in this folder is a synthetic aperture radar (SAR) simulator.
It was written for my thesis called "SAR imaging with a hand-held UWB radar system" at
Blekinge Institute of Technology in Karlskrona, Sweden,
and served the purpose of simulating SAR images, as well as processing real radar images.
Its main purpose was to validate the implementation of common SAR algorithms:

* Chirp waveform generation
* Matched chirp waveform generation
* Pulse-compression
* Radar imaging
* Radio frequency interference generation/simulation
* Radio frequency interference suppression
* SAR imaging (focusing) with global backprojection (GBP)
* SAR 2D FFT generation

It is my hope that the code may be of use to somebody in the future, either as a stepping stone
for continued work, or as something to study to learn how SAR systems behave.

# Code layout
I have tried to keep the code clean, although the simulator has grown out of its boots.
The code is at least segmented into a number of files:

* algorithms.c Contains SAR algorithms.
* file_io.c File input/output (needs to be re-written - it is ugly)
* filters.c Some image filters.
* sar_simulator.c Main program code.
* waveforms.c Generates chirp waveform and its matched signal.
* plot_processed_data.sci Scilab script for plotting processed data (read from file into simulator).
* plot_simulation_data.sci Scilab script for plotting output from simulator run.
* plot_simulation_data.m Matlab script for plotting (contour plots as opposed to 3D plots for the Scilab scripts) output from simulator run.
* sar_simulator.h Structs and function prototypes.
* compile.sh Simple script to compile all source code into the simulator binary "sar_simulator"

# Known bugs
For some frequency and TBP combinations, the simulator crashes in FFTW, and I don't know why.
If anyone finds the reason and can fix it, please e-mail me a patch.

# Required libraries
To compile and link the code you will need the FFTW, math and complex libraries.
The fftw library is located at www.fftw.org
The math and complex libraries are usually included in Linux.
The complex library is needed since the simulator operates on complex numbers.

# Compilation
Before you compile, check the flags set to the C compiler in Makefile
Some of the flags may not be applicable to your processor or compiler version.
The simulator has been successfully compiled using gcc 4.5.3-r2 on Gentoo Linux on an Intel Core i5-2410M CPU.
