CXX = g++
CXXFLAGS = -O3 -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
COMPILE = $(CXX) $(CXXFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: mines

mines: $(OBJFILES)
	$(CXX) -o mines $(OBJFILES) -lgsl -lgslcblas

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	rm -f mines *.o
