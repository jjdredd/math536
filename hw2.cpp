#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "cg.hpp"
#include "util.hpp"

std::vector<double> BVP_Solution(double h, double x1, double x2,
				 double y1, double y2) {

	std::vector<double> res;
	int Nx = (int) floor((x2 - x1)/h);
	int Ny = (int) floor((y2 - y1)/h);
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++){
			double x = x1 + i*h,
				y = y1 + j*h;
			res.push_back(log(x*x + y*y));
		}
	}
	return res;
}

double BVal(int i, int j, double h, double x1, double y1) {
	double x = x1 + i*h,
		y = y1 + j*h;
	return log(x*x + y*y);
}

inline double Space2Mat(unsigned size, int i, int j) {
	return (size - 2) * i + j;
}

inline bool IsBoundX(unsigned N, int i, int j) {
	return (i == 0) || (i == N - 1);
}

inline bool IsBoundY(unsigned N, int i, int j) {
	return (j == 0) || (j == N - 1);
}

// ANYONE ORDERED SPAGHETTI??????
CSM FillMatrix(unsigned N, double h,  double x1, double y1) {

	CSM A(N);
	for (int m = 0; m < N; m++) {
		for (int n = 0; n < N; n++) {
			int row = Space2Mat(N, m, n);
			// y directio fixed, traverse in x direction
			for (int i = m - 1; i <= m + 1; i += 2) {
				if (i < 0 || i > N - 1) continue;
				if (IsBoundX(N, i, n)) continue;
				int col = Space2Mat(N, i, n);
				A.Set();
			}
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				if (j < 0 || j > N - 1) continue;
				if (IsBoundY(N, m, j)) continue;
				int col = Space2Mat(N, m, j);
				A.Set();
			}
			// don't forget the central val

		}
	}
	return A;
}

std::vector<double> SolveLaplaceBVP(double h, double x1, double x2,
				    double y1, double y2, unsigned &steps) {

	int N = (int) floor((x2 - x1)/h);
	if (N != (int) floor((y2 - y1)/h)) {
		std::cout << "Error, Nx != Ny in SolveLaplaceBVP()"
			  << std::endl;
		return res;
	}

	CSM A = FillMatrix();
	std::vector<double> b = FillRHS();

	return CG(A, b, steps, 1e-7);
}

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

	std::cout << "CG: Solving w/o a pc" << std::endl;
	x = CG(A, b, steps, 1e-7);

	std::cout << "ANS:" << std::endl << x << std::endl;

	std::cout << "CG took " << steps << " steps" << std::endl << std::endl;

	std::cout << "CG: Solving w/ a pc = diag(A)" << std::endl;
	x = CGWPC(A, b, steps, 1e-7);

	std::cout << "ANS:" << std::endl << x << std::endl;

	std::cout << "GCWPC took " << steps << " steps" << std::endl;

	return;
}

void Problem_3() {

	return;
}

int main() {

	std::cout << "Problem_2" << std::endl;
	Problem_2();

	return 0;
}
