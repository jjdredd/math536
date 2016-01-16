#include "sbms.hpp"
#include <iostream>


//
// n, m - old cols, rows
// k, l - new cols, rows
// n > m : k = n - m, l = m
// m > n : k = m - n, l = n
// 


SBMS::SBMS(unsigned N, unsigned k) : N(N), k(k) {

	if (k > N) {
		std::cerr << "ERROR: band is wider than the matrix!!"
		     << std::endl;
		return;
	}

	A = new double * [k + 1];
	for (unsigned i = 0; i < k + 1; i++)
		A[i] = new double [N - i];
}

SBMS::~SBMS() {

	for (unsigned i = 0; i < k + 1; i++)
		delete[]  A[i];
	delete[] A;
}

double& SBMS::at(unsigned m, unsigned n) {

	unsigned i, j;
	
	if (m > N || n > N) {
		std::cerr << "OOB in SBMS w/"
			  << "m = " << m
			  << ", n = " << n << std::endl;
		return A[0][0];	// FIXME
	}

	if (n >= m) {
		j = n - m;
		i = m;
	}
	else {
		j = m - n;
		i = n;
	}

	if (j > k || i > N - j) return A[0][0];
	else return A[j][i];
}

double& SBMS::operator() (unsigned m, unsigned n) {
	return at(m, n);
}

SBMS *SBMS::Cholesky() {

	return NULL;
}
