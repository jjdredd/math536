#include <iostream>
#include <fstream>
#include "sbms.hpp"
#include "util.hpp"

// problem 2
void Problem_2() {

	unsigned N = 5;
	SBMS A(N, N);
	std::vector<double> R, b = {5, 3.55, 2.81428571428571,
				    2.34642857142857, 2.01746031746032};
	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = i; j < N; j++)
			A.Set(i, j, 1.0/((i + 1) + (j + 1) - 1));
	}

	std::cout << "Hilbert matrix" << std::endl
		  << A << std::endl;

	LSBMS C = A.Cholesky(); // need smart ptrs

	std::cout << "Cholesky (H)" << std::endl
		  << C << std::endl;


#if 0
	std::cout << "check" << std::endl;

	LSBMS CT = C;
	CT.T();

	std::cout << CT*C << std::endl << std::endl;
#endif

	R = CSolve(C, b);

	std::cout << "x = " << R << std::endl << std::endl;

	std::cout << "perturbed" << std::endl;
	A.Set(0, 4, .20001);

	std::cout << "Perturbed Hilbert matrix" << std::endl
		  << A << std::endl;

	C = A.Cholesky();
	R = CSolve(C, b);

	std::cout << "x = " << R << std::endl;

	return;
}

void Problem_3() {

	// here bandwidth is defined differently
	// than in problem formulation
	unsigned N = 10, k = 3, steps;
	std::vector<double> b = {9, 10, 11,  11, 11, 11, 11,  11, 10, 9};
	SBMS A(N, k);
	ProblemTwoFill(A);
	std::cout << "==========[ A ]==========" << std::endl
		  << A << std::endl << std::endl;

	std::vector<double> x;

	std::ofstream ofile("./conv.txt");
	for (double w = 0.1; w < 2; w += 0.1) {
		x = SORSolve(A, b, steps, w, 1e-6);
		ofile << w << '\t' << steps << std::endl;
	}

	std::cout << "x = " << x << std::endl;

	// std::cout << A * x << std::endl;

}

int main() {

#if 0
	std::cout << "Testing" << std::endl;
	unsigned N = 3;
	SBMS M(N, 2);

	M.Set(2, 2, 0.4);
	
	for (unsigned i = 0; i < N; i++)
		M.Set(i, i, 1.0);

	for (unsigned i = 0; i < N - 1; i++)
		M.Set(i, i + 1, 0.35);

	std::cout << M << std::endl;

	LSBMS *C = M.Cholesky();

	std::cout << *C << std::endl;
#endif

	std::cout << "Problem_2" << std::endl;
	Problem_2();

	std::cout << "Problem_3" << std::endl;
	Problem_3();

	return 0;
}
