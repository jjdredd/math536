#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "util.hpp"
#include "cg.hpp"

inline int Space2Mat(int N, int i, int j) {
	return (N - 2) * i + (j - 1);
}

////////////////////////////////////////////////////////
// Refactor this to eliminate repetitive code !!!!!!! //
////////////////////////////////////////////////////////

CSM FillMatrix(unsigned N, double g) {

	unsigned size = (N - 2) * N;
	CSM A(size);

	for (int m = 1; m < N - 1; m++) {	  // y
		for (int n = 1; n < N - 1; n++) { // x
			int row = Space2Mat(N, m, n);
			// don't forget the central val
			A.Set(row, row, 1 + 2*g);
			// y direction fixed, traverse in x direction
			for (int i = m - 1; i <= m + 1; i += 2) {
				int col = Space2Mat(N, i, n);
				A.Set(row, col, -g/2);
			}
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				if (IsBoundY(N, m, j)) continue;
				int col = Space2Mat(N, m, j);
				A.Set(row, col, -g/2);
			}
		}
	}

	// now take care of boundaries x = 0, x = N - 1
	for (int m = 0; m <= N - 1; m += N - 1) { // x
		for (int n = 1; n < N - 1; n++) { // y
			int col, row = Space2Mat(N, m, n);
			// don't forget the central val
			A.Set(row, row, 1 + 2*g);
			if (!m)	col = Space2Mat(N, m + 1, n);
			else col = Space2Mat(N, m - 1, n);
			A.Set(row, col, -g);
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				if (IsBoundY(N, m, j)) continue;
				int col = Space2Mat(N, m, j);
				A.Set(row, col, -g/2);
			}
		}
	}

	return A;
}

std::vector<double> FillRHS(unsigned N, double g, std::vector<double> v) {

	unsigned size = (N - 2) * N;
	std::vector<double> b(size, 0);

//////////////////////////
// FIXME BOUNDARIES!!!! //
//////////////////////////
// pay attention to y boundary

	for (int m = 1; m < N - 1; m++) {
		for (int n = 1; n < N - 1; n++) {
			int row = Space2Mat(N, m, n);
			int prev = Space2Mat(N, m, n);
			b[row] += (1 - 2*g) * v[prev];
			// y direction fixed, traverse in x direction
			for (int i = m - 1; i <= m + 1; i += 2) {
				prev = Space2Mat(N, i, n);
				b[row] += g/2 * v[prev];
				// we will take care of this bc seperately
			}
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				double prval;
				prev = Space2Mat(N, m, j);
				if (IsBoundY(N, m, j)) prval = j ? 0 : 1;
				else prval = v[prev];
				b[row] += g/2 * prval;
			}
		}
	}

	// now take care of boundaries x = 0, x = N - 1
	for (int m = 0; m <= N - 1; m += N - 1) { // x
		for (int n = 1; n < N - 1; n++) { // y
			int row = Space2Mat(N, m, n);
			int prev = Space2Mat(N, m, n);
			// don't forget the central val
			b[row] += (1 - 2*g) * v[prev];
			// doubled part from ghost nodes
			if (!m)	prev = Space2Mat(N, m + 1, n);
			else prev = Space2Mat(N, m - 1, n);
			b[row] += g * v[prev];
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				double prval;
				prev = Space2Mat(N, m, j);
				if (IsBoundY(N, m, j)) prval = j ? 0 : 1;
				else prval = v[prev];
				b[row] +=  g/2 * prval;
			}
		}
	}

	return b;
}

std::vector<double> StepIBVP(unsigned N, double h, double k,
			     std::vector<double> u, double Error) {
	// u is IC
	CSM A = FillMatrix(N, k / (h * h));
	std::vector<double> b = FillRHS(N, k / (h * h), u);
	unsigned steps;

	return CG(A, b, steps, Error);
}

void PrintSln(unsigned N, double h, std::vector<double> x, std::string file) {

	std::ofstream ofile(file);
	for (unsigned j = 1; j < N - 1; j++) {
		for (unsigned i = 1; i < N - 1; i++) {
			ofile << i * h << '\t' << j * h << '\t'
			      << x[Space2Mat(N, i, j)] << std::endl;
		}
		ofile << std::endl;
	}
}

void Problem_1() {

	std::cout << "Problem_1 (IBVP)" << std::endl;

	unsigned N = 20;
	double Error = 1e-7, e = 1, h = 1.0/N;
	double k = h;
	std::vector<double> u((N - 2) * N, 0);
	std::vector<double> x;
	int n;
	for (n = 0; e > Error; n++) {
		x = StepIBVP(N, h, k, u, Error);
		e = MaxNorm(u - x);
		u = x;		// better use shallow copy or ptr xchg here
	}
	std::cout << "final time " << k * n << ", w/ max error: "
		  << e << " and this took me " << n << " steps" << std::endl;
	PrintSln(N, h, x, "hw4_sln.txt");
}

int main() {

	Problem_1();

	// Problem_2();

	return 0;
}
