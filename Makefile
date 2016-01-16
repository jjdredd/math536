WFLAGS = -W -Wall -Wextra
IFLAGS = -iquote $(shell pwd)/include
OPTFLAGS = -march=native -ffast-math -O3 -funroll-loops # -msse
export CFLAGS = $(IFLAGS) $(WFLAGS) $(OPTFLAGS) -g
export CXXFLAGS = $(CFLAGS) -std=c++11
export LDFLAGS = -lm

SBMS = sbms/sbms.o


hw1.elf: hw1.cpp sbms
	g++ $(CXXFLAGS) $(LDFLAGS) hw1.cpp $(SBMS) -o hw1.elf


sbms:	$(SBMS)
	$(MAKE) -C sbms

clean:
	rm hw1.elf
	$(MAKE) clean -C sbms
