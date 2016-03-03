#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "cg.hpp"
#include "util.hpp"

double BVal(int i, int j, double h, double x1, double y1) {
	double x = x1 + i*h,
		y = y1 + j*h;
	return log(x*x + y*y);
}

inline int Space2Mat(int N, int i, int j) {
	return (N - 2) * (i - 1) + (j - 1);
}

inline bool IsBoundX(int N, int i, int j) {
	return (i == 0) || (i == N - 1);
}

inline bool IsBoundY(int N, int i, int j) {
	return (j == 0) || (j == N - 1);
}

std::vector<double> BVP_Solution(double h, double x1, double x2,
				 double y1, double y2) {

	std::vector<double> res;
	int N = (int) floor((x2 - x1)/h);
	if (N != (int) floor((y2 - y1)/h)) {
		std::cerr << "Error, Nx != Ny in SolveLaplaceBVP()"
			  << std::endl;
	}
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++){
			double x = x1 + i*h,
				y = y1 + j*h;
			res.push_back(log(x*x + y*y));
		}
	}
	return res;
}

CSM FillMatrix(int N) {

	unsigned size = (N - 2) * (N - 2);
	CSM A(size);
	for (int m = 1; m < N - 1; m++) {
		for (int n = 1; n < N - 1; n++) {
			int row = Space2Mat(N, m, n);
			// y directio fixed, traverse in x direction
			for (int i = m - 1; i <= m + 1; i += 2) {
				if (IsBoundX(N, i, n)) continue;
				int col = Space2Mat(N, i, n);
				A.Set(row, col, 1);
			}
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				if (IsBoundY(N, m, j)) continue;
				int col = Space2Mat(N, m, j);
				A.Set(row, col, 1);
			}
			// don't forget the central val
			A.Set(row, row, -4);
		}
	}
	return A;
}

std::vector<double> FillRHS(int N, double h, double x1, double y1) {

	unsigned size = (N - 2) * (N - 2);
	std::vector<double> b(size, 0);
	// suboptimal but easier to write
	for (int m = 1; m < N - 1; m++) {
		for (int n = 1; n < N - 1; n++) {
			int row = Space2Mat(N, m, n);
			// y directio fixed, traverse in x direction
			for (int i = m - 1; i <= m + 1; i += 2) {
				if (!IsBoundX(N, i, n)) continue;
				b[row] -= BVal(i, n, h, x1, y1);
			}
			// x is fixed, traverse in y direction
			for (int j = n - 1; j <= n + 1; j += 2) {
				if (!IsBoundY(N, m, j)) continue;
				b[row] -= BVal(m, j, h, x1, y1);
			}
		}
	}
	return b;
}

std::vector<double> SolveLaplaceBVP(double h, double x1, double x2,
				    double y1, double y2, unsigned &steps) {

	int N = (int) floor((x2 - x1)/h);
	if (N != (int) floor((y2 - y1)/h)) {
		std::cerr << "Error, Nx != Ny in SolveLaplaceBVP()"
			  << std::endl;
	}

	CSM A = FillMatrix(N);
	std::vector<double> b = FillRHS(N, h, x1, y1);

	return CG(A, b, steps, 1e-17);
}

// problem 2
void Problem_2() {

	std::cout << "Problem_2" << std::endl;

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

	char *fname = "./conv_2.txt";
	std::cout << "Problem_3" << std::endl
		  << "output in " << fname << std::endl;

	double h = 0.2;
	double x1 = 1, x2 = 2;
	double y1 = 0, y2 = 1;
	std::ofstream ofile(fname);
	std::vector<double> Nsln;
	for (unsigned i = 0; i < 4; i++, h /= 2) {
		unsigned steps;
		Nsln = SolveLaplaceBVP(h, x1, x2, y1, y2, steps);
		ofile << log(h) << '\t'
		      << log(MaxNorm(Nsln - BVP_Solution(h, x1, x2, y1, y2)))
		      << '\t' << steps << std::endl;
	}

	return;
}

int main() {

	Problem_2();

	Problem_3();

	return 0;
}
