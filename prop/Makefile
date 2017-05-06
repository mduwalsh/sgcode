CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = prop propsim xsum propjob 

.PHONY: all clean

all: $(EXECUTABLES)

prop: prop.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)	

propsim: propsim.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

xsum: xsum.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

propjob: propjob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
