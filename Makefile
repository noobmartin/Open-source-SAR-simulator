CFLAGS=-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -march=native -mcmodel=medium -g
LDFLAGS=-lfftw3 -lm

all: sar_simulator

sar_simulator: waveforms.o algorithms.o sar_simulator.o file_io.o filters.o
	${CC} $^ $(LDFLAGS) -o $@

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@
