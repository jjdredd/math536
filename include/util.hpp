#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <vector>

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	for (unsigned i = 0; i < v.size() - 1; i++) {
		os << v[i] << '\t';
	}
	os << v[v.size() - 1];
	return os;
}

std::vector<double> operator+(std::vector<double>, std::vector<double>);

std::vector<double> operator-(std::vector<double>, std::vector<double>);

double operator*(std::vector<double>, std::vector<double>);

std::vector<double> operator*(std::vector<double>, double);

void PtrXchg(void **, void **);

#endif	// UTIL_HPP
