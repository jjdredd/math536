//
// Banded matrix storage
//



class SBMS {

public:
	SBMS(unsigned, unsigned);
	~SBMS();
	SBMS *Cholesky();
	double& operator()(unsigned, unsigned);

private:
	double& at(unsigned, unsigned);
	unsigned N, k;	   // overall size and band size
	double **A;
};
