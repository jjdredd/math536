#include <cmath>
#include "util.hpp"

std::vector<double> operator+(std::vector<double> lhs,
			      std::vector<double> rhs) {

	std::vector<double> r(rhs.size(), 0);
	if (lhs.size() != rhs.size()) {
		std::cerr << "vector operator+ size mismatch" << std::endl;
		return r;
	}
	for (unsigned i = 0; i < rhs.size(); i++)
		r[i] = lhs[i] + rhs[i];
	return r;
}

std::vector<double> operator-(std::vector<double> lhs,
			      std::vector<double> rhs) {

	std::vector<double> r(rhs.size(), 0);
	if (lhs.size() != rhs.size()) {
		std::cerr << "vector operator+ size mismatch" << std::endl;
		return r;
	}
	for (unsigned i = 0; i < rhs.size(); i++)
		r[i] = lhs[i] - rhs[i];

	return r;
}

double operator*(std::vector<double> lhs, std::vector<double> rhs) {

	double r = 0;
	if (lhs.size() != rhs.size()) {
		std::cerr << "vector operator+ size mismatch" << std::endl;
		return r;
	}
	for (unsigned i = 0; i < rhs.size(); i++)
		r += lhs[i] * rhs[i];
	return r;
}

std::vector<double> operator*(std::vector<double> v, double c) {

	std::vector<double> r(v.size(), 0);
	for (unsigned i = 0; i < v.size(); i++)
		r[i] = v[i] * c;
	return r;
}

std::vector<double> operator*(double c, std::vector<double> v) {

	std::vector<double> r(v.size(), 0);
	for (unsigned i = 0; i < v.size(); i++)
		r[i] = v[i] * c;
	return r;
}

void PtrXchg(void **a, void **b) {
	void *c = *a;
	*a = *b;
	*b = c;
}

// maximum norm used in hw2 problem3
double MaxNorm(std::vector<double> v) {

	double r = 0;
	for (auto a : v) if (fabs(a) > r) r = fabs(a);
	return r;
}
