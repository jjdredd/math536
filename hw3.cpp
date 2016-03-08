#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "util.hpp"

// strictily for prob2
std::vector<double> SolveIBVP(double X, unsigned N, double C, double T) {
	double *u[2];
	std::vector<double> r(N, 0);
	double tau = C * X / N; // v = 1
	int j, n;

	u[0] = new double[N];
	u[1] = new double[N];
	memset(u[0], 0, N * sizeof(double));	// initial values
	u[0][0] = 1;		// left boundary

	for (n = 0; n * tau < T; n++) {
		u[1][0] = 1;	 // left boundary
		u[1][N - 1] = 0; // right boundary
		for (j = 1; j < N - 1; j++) {
			u[1][j] = u[0][j] - C/2 * (u[0][j + 1] - u[0][j - 1])
				+ C*C/2 * (u[0][j + 1] - 2 * u[0][j]
					   + u[0][j - 1]);
		}
		PtrXchg((void **) &u[0], (void **) &u[1]);
	}

	std::copy(u[1], u[1] + N, r.begin());
	delete[] u[0];
	delete[] u[1];

	return r;
}

void Problem_2() {

	std::string base_name("hw3_");
	unsigned N = 100;
	double X = 10, T = 5, h = X / N;
	std::vector<double> Cval = {0.9, 1, 1.1};

	for (double C : Cval){
		std::ofstream of(base_name + std::to_string(C) + ".txt");
		std::vector<double> u = SolveIBVP(X, N, C, T);
		for (unsigned i = 0; i < N; i++)
			of << i * h << '\t' << u[i] << std::endl;
	}
}

int main () {

	Problem_2();

	return 0;
}
