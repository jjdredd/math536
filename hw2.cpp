#include <iostream>
#include <fstream>
#include "cg.hpp"
#include "util.hpp"

// problem 2
void Problem_2() {

	unsigned N = 200, steps;
	CSM A(N);
	std::vector<double> b(N, 0), x;

	for (unsigned i = 0; i < N - 1; i++)
		A.Set(i, i, 1);

	A.Set(N - 1, N - 1, 1e+6);

	for (unsigned i = 0; i < N - 1; i++)
		A.Set(i, i + 1, -0.1);

	for (unsigned i = 0; i < N - 1; i++)
		A.Set(i + 1, i, -0.1);

	b[0] = 0.9;
	b[N - 1] = 999999.9;
	for (unsigned i = 1; i < N - 1; i++) b[i] = 0.8;

	std::cout << CG(A, b, steps, 1e-7) << std::endl;

	return;
}

int main() {

	std::cout << "Problem_2" << std::endl;
	Problem_2();

	return 0;
}
