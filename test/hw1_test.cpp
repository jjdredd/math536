#include <iostream>
#include "sbms.hpp"


int main() {

#if 0
	unsigned N = 3;
	SBMS M(N, 2);


	M.Set(2, 2, 0.4);
	
	for (unsigned i = 0; i < N; i++)
		M.Set(i, i, 1.0);

	for (unsigned i = 0; i < N - 1; i++)
		M.Set(i, i + 1, -0.5);

	std::cout << M;
#endif

	unsigned N = 5;
	SBMS A(N, N);

	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = i; j < N; j++)
			A.Set(i, j, 1.0/((i + 1) + (j + 1) - 1));
	}

	std::cout << "=========[ A ]==========" << std::endl
		  << A << std::endl << std::endl;

	LSBMS C = A.Cholesky();

	LSBMS CT = C;
	CT.T();

	std::cout << "=========[ C ]==========" << std::endl
		  << C << std::endl << std::endl;

	std::cout << static_cast<SBMS> (A*A) << std::endl;


	SBMS B(10, 3);
	ProblemTwoFill(B);
	LSBMS II = (static_cast<LSBMS> (B)).Inv();
	std::cout << II << std::endl << II*B << std::endl;

	return 0;
}
