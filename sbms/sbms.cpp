#include <algorithm>
#include <cmath>
#include "sbms.hpp"


//
// LSBMS Lower Banded Matrix Storage
//

//
// n, m - old cols, rows
// k, l - new cols, rows
// n > m : k = n - m, l = m
// m > n : k = m - n, l = n
// 

void LSBMS::self_alloc() {

	if (k > N) {
		std::cerr << "ERROR: band is wider than the matrix!!"
		     << std::endl;
		return;
	}

	A = new double * [k];
	for (unsigned i = 0; i < k; i++)
		A[i] = new double [N - i];
}

void LSBMS::self_free() {

	for (unsigned i = 0; i < k; i++)
		delete[]  A[i];
	delete[] A;
}

LSBMS::LSBMS(unsigned N, unsigned k) : N(N), k(k), transpose(false) {
	self_alloc();
}

LSBMS::~LSBMS() {
	self_free();
}

double LSBMS::operator() (unsigned m, unsigned n) const {
	return get(m, n);
}

double LSBMS::get(unsigned mm, unsigned nn) const {

	unsigned m = mm, n = nn;
	if (transpose) {
		m = nn;
		n = mm;
	}

	if (m >= N || n >= N) {
		std::cerr << "OOB in LSBMS w/"
			  << "m = " << m
			  << ", n = " << n << std::endl;
		return 0;
	}

	if (n > m) return 0;

	unsigned j = m - n, i = n;

	if (j >= k || i >= N - j) return 0; // out of band
	else return A[j][i];
}

bool LSBMS::Set (unsigned mm, unsigned nn, double a) {

	unsigned m = mm, n = nn;
	if (transpose) {
		m = nn;
		n = mm;
	}

	if (m >= N || n >= N) {
		std::cerr << "Out Of Bounds in LSBMS set("
			  << m << ", " << n << ")" << std::endl;
		return false;
	}

	if (n > m) {
		std::cerr << "Out Of Band in LSMBS set("
			  << m << ", " << n << ")" << std::endl;
		return false;	// upper elements are zero
	}
	
	unsigned j = m - n, i = n;

	if (j >= k || i >= N - j) {
		std::cerr << "Out Of Band in LSMBS set("
			  << m << ", " << n << ")" << std::endl;
		return false;
	}
	else A[j][i] = a;
	return true;
}

void LSBMS::T() {
	transpose = !transpose;
}

LSBMS& LSBMS::operator=(const LSBMS& other){
	if (this == &other) return *this;

	// change the size
	self_free();
	N = other.N;
	k = other.k;
	self_alloc();

	for (unsigned i = 0; i < k; i++)
		std::copy(other.A[i], other.A[i] + N - i, A[i]);
	return *this;
}

std::ostream& operator<<(std::ostream& os, const LSBMS& obj) {

	for (unsigned i = 0; i < obj.N; i++) {
		for (unsigned j = 0; j < obj.N; j++) {
			os << obj.get(i, j) << "\t";
		}
		os << std::endl;
	}
	return os;
}

bool LSBMS::operator==(const LSBMS& other) const {

	std::cerr << " ==  not implemented!" << std::endl;

	for (unsigned i = 0; i < other.N; i++) {
		for (unsigned j = 0; j < other.N; j++) {}
	}
	return true;
}



//
// Symmetric Banded Matrix Storage
//

SBMS::SBMS(unsigned N, unsigned k) : LSBMS(N, k) {}
SBMS::~SBMS() {}

double SBMS::get(unsigned m, unsigned n) const {

	unsigned i, j;

	if (m >= N || n >= N) {
		std::cerr << "OOB in SBMS w/"
			  << "m = " << m
			  << ", n = " << n << std::endl;
		return 0;	// FIXME
	}

	if (n >= m) {
		j = n - m;
		i = m;
	}
	else {
		j = m - n;
		i = n;
	}

	if (j >= k || i >= N - j) return 0;
	else return A[j][i];
}

bool SBMS::Set (unsigned m, unsigned n, double a) {

	unsigned i, j;

	if (m >= N || n >= N) {
		std::cerr << "Out Of Bounds in SBMS set("
			  << m << ", " << n << ")" << std::endl;
		return false;
	}

	if (n >= m) {
		j = n - m;
		i = m;
	}
	else {
		j = m - n;
		i = n;
	}

	if (j >= k || i >= N - j) {
		std::cerr << "Out Of Band in SMBS set("
			  << m << ", " << n << ")" << std::endl;
		return false;
	}
	else A[j][i] = a;
	return true;
}

LSBMS * SBMS::Cholesky() {

	// TODO Optimize for zero elements
	LSBMS *C = new LSBMS(N, N);

	for (unsigned n = 0; n < N; n++) {

		double s = 0;
		for (unsigned i = 0; i < n; i++) 
			s += C->get(n, i) * C->get(n, i);
		
		C->Set(n, n, sqrt(get(n, n) - s));
		
		for (unsigned i = n + 1; i < N; i++) {

			s = 0;
			for (unsigned j = 0; j < n; j++)
				s += C->get(i, j) * C->get(n, j);

			C->Set(i, n, (get(i, n) - s)/C->get(n, n));
		}
		
	}
	
	return C;
}

SBMS operator*(const SBMS A, const SBMS B) {

	// TODO Optimize for zero elements
	unsigned N = A.N;
	SBMS R(N, N);
	if (N != B.N) {
		std::cerr << "MAT MUL SIZE ERROR" << std::endl;
		return R;
	}
	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = 0; j < N; j++) {
			double s = 0;
			for (unsigned n = 0; n < N; n++)
				s += A(i, n) * B(n, j);
			R.Set(i, j, s);
		}
	}
	return R;
}


// CC^Tx = b
// Cy = b
// C^Tx = y
std::vector<double> CSolve(LSBMS &C, std::vector<double>& b) {

	std::vector<double> R(b.size()), y(b.size());
	// Cy = b
	for (int i = 0; i < C.N; i++) {
		double s = 0;
		for (int j = 0; j < i; j++)
			s += C(i, j) * b[j];
		y[i] = (b[i] - s)/C(i, i);
	}
	C.T();
	// C^Tx = y
	for (int i = C.N - 1; i >= 0 ; i--) {
		double s = 0;
		for (int j = C.N - 1; j > i; j--)
			s += C(i, j) * b[j];
		R[i] = (y[i] - s)/C(i, i);
	}
	C.T();
	return R;
}
