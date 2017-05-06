CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = ml allplot fplot mljob xsum sum

.PHONY: all clean

all: $(EXECUTABLES)

ml: ml.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

mljob: mljob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

fplot: fplot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

allplot: allplot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

xsum: xsum.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
	
sum: sum.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
