#include <iostream>
#include "sbms.hpp"

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	for (unsigned i = 0; i < v.size() - 1; i++) {
		os << v[i] << '\t';
	}
	os << v[v.size() - 1];
	return os;
}

// problem 1
void Problem_1() {

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

void Problem_2() {


	// C = L - 1/w * D
	// D = diag(A)
	// U = A - L
	// x_k+1 = -(D/w + L)^(-1)[(1-w^(-1))D + U]x_k + (D/w + L)^(-1)b
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

	std::cout << "Problem_1" << std::endl;
	Problem_1();

	return 0;
}
