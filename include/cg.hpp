#ifndef CG_HPP
#define CG_HPP

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
	CSM(const CSM&);
	~CSM();
	unsigned Size() const;
	void Set(unsigned, unsigned, double);
	double Get(unsigned, unsigned) const;
	CSM& operator=(const CSM&);
	friend std::vector<double> operator* (const CSM,
					      const std::vector<double>);
	friend std::vector<double> operator* (const std::vector<double>,
					      const CSM);

	// obviously, it's better to store row, column numbers and
	// elements in one struct and make an array (vector) of these
	// structs. this will make processor caching much easier.
	// but I'll do as said in the problem and make 2 arrays.
private:
	unsigned N;
	std::vector<double> Elem;
	std::vector<MCoor> C;
};

std::vector<double> CG(const CSM&, std::vector<double>&, unsigned &, double);

std::vector<double> CGWPC(const CSM&, std::vector<double>&, unsigned &, double);

double AProd(const CSM&, std::vector<double>&, std::vector<double>&);

double ANorm(const CSM&, std::vector<double>&);

#endif	// CG_HPP
