CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = mlp xsum imghtml

.PHONY: all clean

all: $(EXECUTABLES)

mlp: mlp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
	
xsum: xsum.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
	
imghtml: imghtml.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
