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

std::vector<double> operator*(double, std::vector<double>);

double MaxNorm(std::vector<double>);

void PtrXchg(void **, void **);

inline int Space2Mat(int N, int i, int j) {
	return (N - 2) * (i - 1) + (j - 1);
}

inline bool IsBoundX(int N, int i, int j) {
	return (i == 0) || (i == N - 1);
}

inline bool IsBoundY(int N, int i, int j) {
	return (j == 0) || (j == N - 1);
}

#endif	// UTIL_HPP
