CC = gcc
CFLAGS = -Wall -fopenmp -std=c99
OBJS = initialize.o density.o io.o main.o memory.o vars.o powerspectrum.o 
INCLS = initialize.h density.h io.h	memory.h powerspectrum.h vars.h
EXEC = POWERSPEC
FFTW_INCL = -I/home/rui/library/fftw3/include
FFTW_LIBS = -L/home/rui/library/fftw3/lib

$(EXEC) : $(OBJS)
	$(CC) $(CFLAGS) $(FFTW_LIBS) $(OBJS) -lfftw -lm -o $(EXEC)

$(OBJS) : $(INCLS)
	$(CC) $(CFLAGS) $(FFTW_INCL) -c -o $*.o $*.c

clean :
	rm -f $(EXEC) $(OBJS)