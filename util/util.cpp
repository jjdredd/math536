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
