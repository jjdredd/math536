#include <iostream>
#include <vector>

//
// LSBMS Lower Symmetric Banded Matrix Storage
//

class LSBMS {

public:
	LSBMS(unsigned, unsigned);
	LSBMS(const LSBMS&);
	~LSBMS();
	double operator()(unsigned, unsigned) const;
	LSBMS& operator=(const LSBMS&);
	bool operator==(const LSBMS&) const;
	LSBMS operator-(const LSBMS&) const;
	LSBMS operator+(const LSBMS&) const;
	LSBMS Inv();
	virtual bool Set(unsigned, unsigned, double);
	virtual double get(unsigned, unsigned) const;
	void T();
	friend std::ostream& operator<<(std::ostream&, const LSBMS&);
	std::vector<double> Solve(std::vector<double>&);
	LSBMS operator*(const LSBMS&) const;
	LSBMS operator*(const double) const;
	std::vector<double> operator*(const std::vector<double>&) const;

protected:
	void self_alloc();
	void self_free();
	unsigned N, k;	   // overall size and band size
	double **A;
	bool transpose;
};


//
// Symmetric Banded Matrix Storage
//

class SBMS : public LSBMS {

public:
	SBMS(unsigned, unsigned);
	SBMS(const LSBMS&);
	SBMS(const SBMS&);
	~SBMS();
	LSBMS Cholesky();
	virtual bool Set(unsigned, unsigned, double);
	virtual double get(unsigned, unsigned) const;
	friend std::vector<double> SORSolve(SBMS&, std::vector<double>&,
					    unsigned *,
					    double w = 1, double e = 1e-6);
	friend void ProblemTwoFill(SBMS &);
};


std::vector<double> CSolve(LSBMS&, std::vector<double>&);
