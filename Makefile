CFLAGS=-O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -march=native -mcmodel=medium -g
LDFLAGS=-lfftw3 -lm

all: sar_simulator

sar_simulator: waveforms.o file_io.o algorithms.o sar_simulator.o filters.o common_functions.o
	${CC} $(wildcard obj/*.o) $(LDFLAGS) -o $@

%.o: **/%.c
	${CC} ${CFLAGS} -c $< -o obj/$@

clean:
	rm obj/*
	rm sar_simulator
