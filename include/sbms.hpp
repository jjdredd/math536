#include <iostream>
//
// LSBMS Lower Symmetric Banded Matrix Storage
//

class LSBMS {

public:
	LSBMS(unsigned, unsigned);
	~LSBMS();
	double operator()(unsigned, unsigned) const;
	LSBMS& operator=(const LSBMS&);
	virtual bool Set(unsigned, unsigned, double);
	virtual double get(unsigned, unsigned) const;
	friend std::ostream& operator<<(std::ostream&, const LSBMS&);

protected:
	void self_alloc();
	void self_free();
	unsigned N, k;	   // overall size and band size
	double **A;
};


//
// Symmetric Banded Matrix Storage
//

class SBMS : public LSBMS {

public:
	SBMS(unsigned, unsigned);
	~SBMS();
	LSBMS *Cholesky();
	virtual bool Set(unsigned, unsigned, double);
	virtual double get(unsigned, unsigned) const;
};
