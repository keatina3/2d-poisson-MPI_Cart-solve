CC = mpicc

CFLAGS = -g -Wall

LDFLAGS = -lm 

OBJECTS = jacobi.o prog.o poisson.o utils.o

prog: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test: prog
	mpiexec -n 16 ./prog -g 9

.PHONY: clean test

clean:
	$(RM) *.o prog
