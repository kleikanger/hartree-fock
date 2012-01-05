
class HFmatrix{

	//Member variables
	double** ppHFdata;
	double**** coulombIntegrals;
	int nCutoff;
   	int iNumberOfParticles; //usually=2nCutoff
	int iNrMeshpt;
    double dIntLimMin;
    double dIntLimMax;
	double* pOrbitalEnergies;;
	double** ddUnitaryMatrix;
	
	//Methods

	public:
	//member variables

	//Methods (member functions)
	HFmatrix(int);
	void transform();
	void diagonalize(double*);
	double findEnergy();	

};

