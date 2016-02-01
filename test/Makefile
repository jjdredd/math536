
SBMS = ../sbms/sbms.o
UTIL = ../util/util.o

all: hw1_test.elf

hw1_test.elf: hw1_test.cpp $(SBMS) $(UTIL)
	g++ $(CXXFLAGS) $(LDFLAGS) hw1_test.cpp $(UTIL) $(SBMS) -o hw1_test.elf

clean:
	rm hw1_test.elf
