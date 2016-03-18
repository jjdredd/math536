#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "util.hpp"
#include "cg.hpp"


std::vector<double> StepIBVP(unsigned N, double k, std::vector<double> u
			     double Error) {
	CSM A = FillMatrix(N);
	std::vector<double> b = FillRHS(N, h, x1, y1);

	return CG(A, b, steps, Error);
}

void Problem_1() {

	unsigned N = 20;
	double Error = 1e-17, e = 1, k = 1.0/N;
	std::vector<double> u(N, 0);
	for (int n = 0; e > Error; n++) {
		std::vector<double> x = StepIBVP(N, k, u, Error);
		e = MaxNorm(u - x);
		u = x;		// better use shallow copy or ptr exchange here
	}
}

int main() {

	Problem_1();

	// Problem_2();

	return 0;
}
