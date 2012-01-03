
class HFmatrix{

	//Member variables
	double** ppHFdata;
	double**** coulombIntegrals;
	int nCutoff;
   	int iNumberOfParticles; //usually=2nCutoff
	int iNrMeshpt;
    double dIntLimMin;
    double dIntLimMax;

	public:
	//member variables
	double* pOrbitalEnergies;

	//Methods (member functions)
	HFmatrix(int);
	void transform(double**);
	void diagonalize(double*, double**);
	
};

