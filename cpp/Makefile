CXX = g++
CXXFLAGS = -O3
COMPILE = $(CXX) $(CXXFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: mines

mines: $(OBJFILES)
	$(CXX) -o mines $(OBJFILES) -l gsl -l gslcblas

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	rm -f mines *.o
