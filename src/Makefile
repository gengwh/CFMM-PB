#
# linux
LIBSBLAS = -lblas -llapack
LIBMATH = -lm
LIBF2C = -lf2c

CFLAGS = -O2
#CFLAGS  = -g
CC = gcc -c
LD = gcc


HEADERS = gkGlobal.h gk.h f2c.h
OBJ = gkSetup.o input.o gkGlobal.o expan.o \
      moments.o numQuad.o fmm.o kernel.o \
      precond_fmm.o gmres.o treecode.o

%.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) -o $@ $<

coulomb: coulomb.o $(OBJ)
	$(LD) -o coulomb coulomb.o $(OBJ) $(LIBMATH) $(LIBSBLAS)

clean:
	\rm -f *.o coulomb
