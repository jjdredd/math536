
SBMS = ../sbms/sbms.o

all: hw1_test.elf

hw1_test.elf: hw1_test.cpp $(SBMS)
	g++ $(CXXFLAGS) $(LDFLAGS) hw1_test.cpp $(SBMS) -o hw1_test.elf

clean:
	rm hw1_test.elf
