#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "util.hpp"
#include "cg.hpp"


//////////////////////////////
// THIS CODE IS A BIT MESSY //
//////////////////////////////


double ExactSln(double x) {
	return (1 / (1 - x/2) - 1) / 2;
}

std::vector<double> FillESln(unsigned N, double h) {

	std::vector<double> b(N - 1, 0);

	for (unsigned i = 1; i < N; i++) b[i - 1] = ExactSln(i * h);

	return b;
}

double Error(const std::vector<double> &v) {

	double e = 0;
	for (auto &x : v) e += x * x;

	return sqrt(e/v.size());
}

double A(double x) {
	return (1 - x/2) * (1 - x/2);
}

// l'
double dl(unsigned n, double h, double x) {
	if (x > (n + 1) * h || x < (n - 1) * h) return 0;
	if (x > n * h) return -1/h;
	return 1/h;
}

// map [-1,1] into [a,b]
double coor_map(double a, double b, double x) {
	if (b <= a) std::cerr << "ERROR: b <= a" << std::endl;
	return ((b - a) * x + a + b) / 2;
}

// second gauss quad
// careful! i is the 'central' index
double GaussQuad(unsigned N, unsigned i, unsigned j, double h) {

	double b = std::min(i + 1, N - 1) * h;
	double a = (i - 1) * h;
	double x1 = coor_map(a, b, -1/sqrt(3));
	double x2 = coor_map(a, b, 1/sqrt(3));

	return (b - a)/2 * (A(x1) * dl(i, h, x1) * dl(j, h, x1)
			    + A(x2) * dl(i, h, x2) * dl(j, h, x2));
}

// fill galerkin matrix
CSM FillGMat(unsigned N, double h) {

	CSM A(N - 1);
	for (unsigned i = 1; i < N; i++) {
		for (unsigned j = i - 1; j <= i + 1 ; j++) {
			// i, j are coordinates on grid
			// need to transform (shift) them to
			// coordinates in the matrix
			if (j > 0 && j < N) {
				A.Set(i - 1, j - 1, GaussQuad(N, i, j, h));
			}
		}
	}
	return A;
}

void Problem_3() {

	std::string file("hw5_errors.txt");
	std::cout << "Problem_3, ouput is in " << file << std::endl;

	std::ofstream ofile(file);
	double e = 1e-7;
	unsigned steps;

	for (double h = 0.1; h >= 0.025; h /= 2) {
		unsigned N = static_cast<unsigned> (floor(1/h));
		// rhs
		std::vector<double> b(N - 1, 0);
		b[N - 2] = A(N * h); // right bc
		CSM M = FillGMat(N, h);
		std::vector<double> x = CG(M, b, steps, e);
		ofile << log(h) << '\t' << log(Error(x - FillESln(N, h)))
		      << std::endl;
	}
}

int main() {

	Problem_3();

	return 0;
}
