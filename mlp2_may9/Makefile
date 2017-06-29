CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = mlp mlpjob 

.PHONY: all clean

all: $(EXECUTABLES)

mlp: mlp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
	
mlpjob: mlpjob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
