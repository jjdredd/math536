#include <iostream>
#include "sbms.hpp"


int main() {

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

	return 0;
}
