WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops # -msse
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS) -g
export CXXFLAGS = $(CFLAGS) -std=c++11
export LDFLAGS = -lm

SBMS = sbms/sbms.o
UTIL = util/util.o
CG = cg/cg.o

all: hw1.elf hw2.elf hw3.elf test

hw1.elf: hw1.cpp sbmslib utillib
	g++ $(CXXFLAGS) $(LDFLAGS) hw1.cpp $(SBMS) $(UTIL) -o hw1.elf

hw2.elf: hw2.cpp cglib utillib
	g++ $(CXXFLAGS) $(LDFLAGS) hw2.cpp $(CG) $(UTIL) -o hw2.elf

hw3.elf: hw3.cpp utillib
	g++ $(CXXFLAGS) $(LDFLAGS) hw3.cpp $(UTIL) -o hw3.elf

hw4.elf: hw4.cpp utillib cglib
	g++ $(CXXFLAGS) $(LDFLAGS) hw4.cpp $(CG) $(UTIL) -o hw4.elf

test:	sbmslib utillib
	$(MAKE) -C test

sbmslib:
	$(MAKE) -C sbms

utillib:
	$(MAKE) -C util

cglib: utillib
	$(MAKE) -C cg

clean:
	-$(MAKE) clean -C sbms
	-$(MAKE) clean -C test
	-$(MAKE) clean -C util
	-rm hw1.elf
	-rm hw2.elf
	-rm hw3.elf
