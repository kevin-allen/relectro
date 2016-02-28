# Makefile for building the C language shared library for the Relectro
C = gcc
OPTS = -c -fPIC
LOADER = gcc

OBJECTS = autocorrelation.c

 CalcMatrixAvg.so: $(OBJECTS)
	R CMD SHLIB -o relectro.so $(OBJECTS)

.c.o: ; $(C) $(OPTS) -c $<

clean:
	-rm *.o *.so
