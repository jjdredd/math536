#include "cg.hpp"
#include "util.hpp"
#include <cmath>

CSM::CSM(unsigned size) : N(size) {}

CSM::CSM(const CSM& other) {
	N = other.N;
	C = other.C;
	Elem = other.Elem;
}

CSM::~CSM() {}

unsigned CSM::Size() const {
	return N;
}

void CSM::Set(unsigned row, unsigned col, double val) {

	if (row >= N || col >= N) {
		std::cerr << "row or column out of range ("
			  << row << ", " << col << ")"
			  << std::endl;
		return;
	}
	MCoor mc = {row, col};
	Elem.push_back(val);
	C.push_back(mc);
}

double CSM::Get(unsigned m, unsigned n) const {

	for (unsigned ind = 0; ind < C.size(); ind++) {
		if (m == C[ind].row && n == C[ind].col)	 return Elem[ind];
	}
	return 0;
}

CSM& CSM::operator=(const CSM& other) {

	if (this == &other) return *this;
	N = other.N;
	C = other.C;
	Elem = other.Elem;

	return *this;
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

std::vector<double> operator* (const std::vector<double> b,
				      const CSM A) {

	std::cerr << "ERROR: vector x CSM not implemented!!!" << std::endl;
	// FIXME: rewrite this function. Needs a fix.
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

std::ostream& operator<<(std::ostream& os, const CSM& A) {
	for (unsigned i = 0; i < A.N; i++) {
		for (unsigned j = 0; j < A.N; j++) {
			os << A.Get(i, j) << '\t';
		}
		os << std::endl;
	}
	return os;
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
		error = sqrt(r * r) / sqrt(N);
		steps++;
	} while (error > e);

	return u;
}

std::vector<double> CGWPC(const CSM& A, std::vector<double>& b,
		       unsigned &steps, double e) {

	unsigned N = b.size();
	double a, error;
	std::vector<double> r, p, u(N, 0);
	CSM BI(N);

	for (unsigned i = 0; i < N; i++) {
		BI.Set(i, i, 1 / A.Get(i, i));
	}

	steps = 0;
	r = b - A * u;

	p = BI * r;

	do {
		a = (1 / ANorm(A, p)) * r * p;
		u = u + a * p;
		r = r - a * (A * p);
		std::vector<double> z = BI * r;
		p = z - (1 / ANorm(A, p)) * AProd(A, z, p) * p;
		error = sqrt(r * r);
		steps++;
	} while (error > e);

	return u;
}
