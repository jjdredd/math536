#include <vector>
#include <iostream>

//
// compact storage matrix class
//

struct MCoor {
	unsigned row, col;
};

class CSM {

public:
	CSM(unsigned);
	~CSM();
	Set(unsigned, unsigned, double);
	friend CSM operator* (const CSM, const std::vector<double>&);

	// obviously, it's better to store row, column numbers and
	// elements in one struct and make an array (vector) of these
	// structs. this will make processor caching much easier.
	// but I'll do as said in the problem and make 2 arrays.
private:
	std::vector<double> Elem;
	std::vector<MCoor> C;
	unsigned N;
};

std::vector<double> CG(const CSM&, std::vector<double>& b, int &steps);
