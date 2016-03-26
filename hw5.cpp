#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "util.hpp"
#include "cg.hpp"

double ExactSln(double x) {
	return 1/2 * (1 / (1 - x/2) - 1);
}

double A(double x) {
	return (1 - x/2) * (1 - x/2);
}

double dl(unsigned n, double h, double x) {
	if (x > (n + 1) * h || x < (n - 1) * h) return 0;
	if (x > n * h) return -1;
	else return 1;
}

// map [a,b] into [-1,1]
double coor_map(double a, double b, double x) {
	if (b <= a) std::cerr << "ERROR: b <= a" << std::endl;
	return ((b - a) * x + a + b) / 2;
}

double GaussQuad(unsigned i, unsigned j, double h) {
	double x1 = coor_map(-1/sqrt(3));
	double x2 = coor_map(1/sqrt(3));
		return () * A();
}
