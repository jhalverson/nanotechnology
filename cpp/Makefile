CC = g++
CFLAGS = -O3
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: mines

mines: $(OBJFILES)
	$(CC) -o mines $(OBJFILES) -l gsl -l gslcblas

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	rm -f mines *.o
