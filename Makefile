WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops # -msse
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS) -g
export CXXFLAGS = $(CFLAGS) -std=c++11
export LDFLAGS = -lm

SBMS = sbms/sbms.o


hw1.elf: hw1.cpp sbmslib
	g++ $(CXXFLAGS) $(LDFLAGS) hw1.cpp $(SBMS) -o hw1.elf


sbmslib:
	$(MAKE) -C sbms

test:	sbmslib
	$(MAKE) -C test

clean:
	rm hw1.elf
	$(MAKE) clean -C sbms
	$(MAKE) clean -C test
