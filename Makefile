CC = mpicc

CFLAGS = -g -Wall

LDFLAGS = -lm 

OBJECTS = jacobi.o prog.o poisson.o utils.o

prog: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) *.o prog
