CC = mpicc

CFLAGS = -g -Wall

OBJECTS = 

prog: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) *.o prog
