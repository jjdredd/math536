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
	for (unsigned i = 0; i < k; i++) {
		A[i] = new double [N - i];
		std::fill_n(A[i], N - i, 0);
	}

}

void LSBMS::self_free() {

	for (unsigned i = 0; i < k; i++)
		delete[]  A[i];
	delete[] A;
}

LSBMS::LSBMS(unsigned N, unsigned k) : N(N), k(k), A(NULL), transpose(false) {
	self_alloc();
}

LSBMS::LSBMS(const LSBMS& other) {
	// change the size
	transpose = other.transpose;
	N = other.N;
	k = other.k;
	self_alloc();

	for (unsigned i = 0; i < k; i++)
		std::copy(other.A[i], other.A[i] + N - i, A[i]);
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
	transpose = other.transpose;
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

LSBMS LSBMS::operator*(const LSBMS& B) const {

	if (transpose != B.transpose && std::min(k, B.k) != 1) {
		std::cerr << "don't multiply L and U" << std::endl;
	}
	// TODO Optimize for zero elements
	LSBMS R(N, N);
	if (N != B.N) {
		std::cerr << "MAT MUL SIZE ERROR" << std::endl;
		return R;
	}
	if (!transpose) {
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j <= i ; j++) {
				double s = 0;
				for (unsigned n = 0; n < N; n++)
					s += get(i, n) * B(n, j);
				R.Set(i, j, s);
			}
		}
	} else {
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = i; j < N ; j++) {
				double s = 0;
				for (unsigned n = 0; n < N; n++)
					s += get(i, n) * B(n, j);
				R.Set(i, j, s);
			}
		}
	}

	return R;
}

LSBMS LSBMS::operator*(const double a) const {

	LSBMS R(N, k);
	for (unsigned i = 0; i < k; i ++) {
		for (unsigned j = 0; j < N - i; j++) {
			R.A[i][j] = a * A[i][j];
		}
	}
	return R;
}

// can be optimized
std::vector<double> LSBMS::operator*(const std::vector<double>& b) const {

	if (b.size() != N) std::cerr << "ERROR: wrong vector/matrix size"
				     << std::endl;

	std::vector<double> r(N, 0);
	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = 0; j < N; j++) {
			r[i] += get(i, j) * b[j];
		}
	}
	return r;
}

LSBMS LSBMS::operator-(const LSBMS& o) const {

	unsigned kk = std::min(k, o.k);
	LSBMS R(N, k);

	if (k != kk) R = *this;
	else R = o;

	if (transpose != o.transpose && kk != 1)
		std::cerr << "operator- can't use L and U" << std::endl;

	if (N != o.N)
		std::cerr << "operator- size mismatch" << std::endl;

	for (unsigned i = 0; i < kk; i ++) {
		for (unsigned j = 0; j < N - i; j++) {
			R.A[i][j] = A[i][j] - o.A[i][j];
		}
	}

	return R;
}

LSBMS LSBMS::operator+(const LSBMS& o) const {

	unsigned kk = std::min(k, o.k);
	LSBMS R(N, k);

	if (k != kk) R = *this;
	else R = o;

	if (transpose != o.transpose && kk != 1)
		std::cerr << "operator+ can't use L and U" << std::endl;

	if (N != o.N)
		std::cerr << "operator+ size mismatch" << std::endl;

	for (unsigned i = 0; i < kk; i ++) {
		for (unsigned j = 0; j < N - i; j++) {
			R.A[i][j] = A[i][j] + o.A[i][j];
		}
	}

	return R;

}

std::vector<double> LSBMS::Solve(std::vector<double>& b) {

	std::vector<double> R(b.size());
	if (!transpose) {
		// Cy = b
		for (int i = 0; i < N; i++) {
			double s = 0;
			for (int j = 0; j < i; j++)
				s += get(i, j) * R[j];
			R[i] = (b[i] - s)/get(i, i);
		}
	} else {	// transpose
		// C^Tx = y
		for (int i = N - 1; i >= 0 ; i--) {
			double s = 0;
			for (int j = N - 1; j > i; j--)
				s += get(i, j) * R[j];
			R[i] = (b[i] - s)/get(i, i);
		}
	}
	return R;
}

LSBMS LSBMS::Inv() {

	LSBMS Res(N, N);
	for (unsigned n = 0; n < N; n++) {
		std::vector<double> I(N, 0), r;
		I[n] = 1;
		r = Solve(I);
		for (unsigned j = 0; j < N; j++)
			Res.Set(j, n, r[j]);
       }

	return Res;
}

//
// Symmetric Banded Matrix Storage
//

SBMS::SBMS(unsigned N, unsigned k) : LSBMS(N, k) {}

SBMS::SBMS(const LSBMS& o) : LSBMS(o) {}

SBMS::SBMS(const SBMS& o)  : LSBMS(o) {}

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

LSBMS SBMS::Cholesky() {

	// TODO Optimize for zero elements
	LSBMS C(N, N);

	for (unsigned n = 0; n < N; n++) {

		double s = 0;
		for (unsigned i = 0; i < n; i++) 
			s += C.get(n, i) * C.get(n, i);
		
		C.Set(n, n, sqrt(get(n, n) - s));
		
		for (unsigned i = n + 1; i < N; i++) {

			s = 0;
			for (unsigned j = 0; j < n; j++)
				s += C.get(i, j) * C.get(n, j);

			C.Set(i, n, (get(i, n) - s)/C.get(n, n));
		}
		
	}
	
	return C;
}


// CC^Tx = b
// Cy = b
// C^Tx = y
std::vector<double> CSolve(LSBMS &C, std::vector<double>& b) {

	std::vector<double> R, y;

	y = C.Solve(b);
	C.T();
	R = C.Solve(y);
	C.T();

	return R;
}

//
// C = L - 1/w * D
// D = diag(A)
// U = A - L
// x_k+1 = -(D/w + L)^(-1)[(1-w^(-1))D + U]x_k + (D/w + L)^(-1)b
//
std::vector<double> SORSolve(SBMS& A, std::vector<double>& b, unsigned *steps,
			     double w, double e) {

	double error = 0;
	unsigned N = A.N, k = A.k;
	std::vector<double> r(N, 0),
		*x1 = new std::vector<double>,
		*x2 = new std::vector<double>;
	LSBMS D(N, 1), L(N, k), U(N, k);

	x1->assign(N, 0);
	x2->assign(N, 0);
	for (unsigned i = 0; i < N; i++) D.Set(i, i, A.get(i, i));
	for (unsigned i = 0; i < N; i++) {
		for (unsigned j = 0; j < i; j++) {
			L.Set(i, j, A.get(i, j));
		}
	}
	U = L;
	U.T();

	// std::cout << std::endl
	// 	  << "============[ D ]===========" << std::endl
	// 	  << D << "============[ L ] ============" << std::endl
	// 	  << L << "============[ U ] ============" << std::endl
	// 	  << U << std::endl << std::endl;

	LSBMS Ml = (D * (1/w) + L).Inv(), Mu = D * (1 - 1/w) + U;

	// std::cout << std::endl
	// 	  << "============[ Ml ]===========" << std::endl
	// 	  << Ml << "============[ Mu ] ============" << std::endl
	// 	  << Mu << std::endl;

	*steps = 0;
	do {
		std::vector<double> a = Mu * (*x1), *auxptr, diff(N, 0);

		for (unsigned i = 0; i < N; i++)
			diff[i] = b[i] - a[i];
		*x2 = Ml * diff;
		auxptr = x2;
		x2 = x1;
		x1 = auxptr;
		error = 0;
		for (unsigned i = 0; i < N; i++) {
			error += (x1->at(i) - x2->at(i))
				* (x1->at(i) - x2->at(i));
		}
		error = sqrt(error);
		// std::cout << error << std::endl;
		(*steps)++;
	} while (error > e);

	std::copy(x1->begin(), x1->end(), r.begin());
	delete x1;
	delete x2;

	return r;
}

void ProblemTwoFill(SBMS &A) {

	std::fill_n(A.A[0], A.N, 7);
	std::fill_n(A.A[1], A.N - 1, 1);
	std::fill_n(A.A[2], A.N - 2, 1);
}
