#include "cg.hpp"

CSM::CSM(unsigned size) : N(size) {}

CSM::~CSM() {}

void CSM::Set(unsigned row, unsigned col, double val) {

	MCoor mc = {row, col};
	Elem.push_back(val);
	C.push_back(mc);
}

std::vector<double> operator* (const CSM, const std::vector<double>) {

	std::vector<double> r;

	return r;
}
