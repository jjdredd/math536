WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops # -msse
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS) -g
export CXXFLAGS = $(CFLAGS) -std=c++11
export LDFLAGS = -lm

SBMS = sbms/sbms.o
UTIL = util/util.o

all: hw1.elf test

hw1.elf: hw1.cpp sbmslib utillib
	g++ $(CXXFLAGS) $(LDFLAGS) hw1.cpp $(SBMS) $(UTIL) -o hw1.elf

test:	sbmslib utillib
	$(MAKE) -C test

sbmslib:
	$(MAKE) -C sbms

utillib:
	$(MAKE) -C util

clean:
	-$(MAKE) clean -C sbms
	-$(MAKE) clean -C test
	-$(MAKE) clean -C util
	-rm hw1.elf
