# Contact
You may contact me (Alexander Rajula) at alexander@rajula.se or superrajula@gmail.com

If you have any questions about the code, or have found it useful, don't hesitate to e-mail me.

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
