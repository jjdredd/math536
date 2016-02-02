#include "cg.hpp"
#include "util.hpp"
#include <cmath>

CSM::CSM(unsigned size) : N(size) {}

CSM::~CSM() {}

unsigned CSM::Size() const {
	return N;
}

void CSM::Set(unsigned row, unsigned col, double val) {

	if (row >= N || col >= N) {
		std::cerr << "row or column out of range" << std::endl;
		return;
	}
	MCoor mc = {row, col};
	Elem.push_back(val);
	C.push_back(mc);
}

std::vector<double> operator* (const CSM A, const std::vector<double> b) {

	if (b.size() != A.N) {
		std::cerr << "size mismatch in operator*" << '\t'
			  << b.size() << '\t' << A.N << std::endl;
	}
	std::vector<double> r(b.size(), 0);
	for (unsigned i = 0; i < A.N; i++) {
		for (unsigned ind = 0; ind < A.C.size(); ind++) {
			if (A.C[ind].row == i) {
				r[i] += A.Elem[ind] * b[ A.C[ind].col ];
			}
		}
	}
	return r;
}

double AProd(const CSM& A, std::vector<double>& a, std::vector<double>& b) {
	return a * (A * b);
}

double ANorm(const CSM& A, std::vector<double>& a) {
	return a * (A * a);
}

std::vector<double> CG(const CSM& A, std::vector<double>& b,
		       unsigned &steps, double e) {

	unsigned N = b.size();
	double a, error;
	std::vector<double> r, p, u(N, 0);

	steps = 0;
	p = r = b - A * u;

	do {
		a = (1 / ANorm(A, p)) * r * p;
		u = u + a * p;
		r = r - a * (A * p);
		p = r - (1 / ANorm(A, p)) * AProd(A, r, p) * p;
		error = sqrt(r * r);
		steps++;
	} while (error > e);

	return u;
}
