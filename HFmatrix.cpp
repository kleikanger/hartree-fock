/*

Some parameters defined in the constructor HFmatrix::HFmatrix. Should be taken as argument.
Class methods, objects, .... called from main. Can be moved to other program.

   */

//if PRINT == true some results will be printed
#define PRINT true

#include <iostream>
#include "lib/lib.h"
#include "HFmatrix.h"
#include <cmath>

double calculateInnerIntegrals(double , int , int );  //REMOVE!

using std::cout;

//Function not programmed by me (MH-J?).
double LaguerreGeneral(int,double,double);
//Functions programmed by me
double radialWF(double,int,int);
double hFMatrixElements(int ,int ,int ,int ,int ,int ,double, double);

//Constructs Hartre-Fock (iNumberOfParticles x iNumberOfParticles ) martix
HFmatrix::HFmatrix(int nCutoff){ //XXX May need more parameters as input. 
//startvimfold

	nCutoff=2; // XXX remove
	
	//XXX Should may be part of constructor call?
	iNumberOfParticles=2*nCutoff;
	//There might be a lot of round off error with to many terms because of roundofferror in the integration routine.	
	//For some reason, david's integrationroutine with two separate inner and outer integrals converged much faster (?).
	
	//Calc of matrix elements O(iNrMeshpt**2) 
	iNrMeshpt=1500;
	dIntLimMin=0.0;
	dIntLimMax=250.0;
	
	double **ddUnity;

    //Declare elements of variable:coulombIntegrals:
    coulombIntegrals = new double***[nCutoff];
	    for (int i = 0; i < nCutoff; i++) {
	        coulombIntegrals[i] = new double**[nCutoff];
	        for (int j = 0; j < nCutoff; j++) {
	            coulombIntegrals[i][j] = new double*[nCutoff];
	            for (int k = 0; k < nCutoff; k++) {
	                coulombIntegrals[i][j][k] = new double[nCutoff];
	                for (int l = 0; l < nCutoff; l++) {
	                    coulombIntegrals[i][j][k][l] = 0;
	                }
		 		}
			}
		}	
	
	/*Initialize elements of variable:coulombIntegrals
  	finds hartree and fock part of the  matrixelement (<ab|cd>) by integratian of the direct and exchange term.
  	remember that interchange of the terms a,c and b,d will give the same results. therefore only some of the 
	elements of the  matrix will be neccasary to compute. this fact can be shown analytically (easily).
 	*/
	    for (int n1 = 0; n1 < nCutoff; n1++) {
			for(int n2 = 0; n2 < nCutoff; n2++) {
	            for(int n3 = n1; n3 < nCutoff; n3++) {
	                for(int n4 = n2; n4 < nCutoff; n4++) {
                    coulombIntegrals[n3][n2][n1][n4] =
                    coulombIntegrals[n1][n4][n3][n2] =
                    coulombIntegrals[n3][n4][n1][n2] =
                    coulombIntegrals[n1][n2][n3][n4] =
					hFMatrixElements(n1+1, n2+1, n3+1, n4+1,iNumberOfParticles,iNrMeshpt,dIntLimMin,dIntLimMax);
					cout<<coulombIntegrals[n3][n2][n1][n4]<<"\t";
					}
					cout<<".\n";
			    }
		 	}
		}
	
	//Declaring and constructing unity matrix	
  	ddUnity = (double **) matrix( iNumberOfParticles, iNumberOfParticles , sizeof(double));
	for (int a=0; a< iNumberOfParticles; a++){
		ddUnity[a][a]=1.0;
		for (int b=a+1; b< iNumberOfParticles; b++){
 			ddUnity[a][b]=ddUnity[b][a]=0.0;
		}
	}
	//initializing eigenvalues(energies) of orbitals
	pOrbitalEnergies = new double[iNumberOfParticles];
	for (int i = 0; i < iNumberOfParticles/2; i++) {
		pOrbitalEnergies[i] = -((double) (iNumberOfParticles*iNumberOfParticles))/((double) (2*(i+1)*(i+1))); //iNumberOfParticles really Z!!
	}
    
	//Declare elements of variable:ppHFdata
 	ppHFdata= (double **) matrix( iNumberOfParticles, iNumberOfParticles, sizeof(double));
	
	//Finding initial matrix
	transform(ddUnity);
	
	#if A
	//Prin test
	for (int i=0;i<4;i++){
		for (int j=0;j<4;j++){
		cout<<", \t"<<ppHFdata[i][j];
		}
	cout<<"\n";
	}
	#endif
	//XXX SLETT ARRAYS!

}//END of constructor HFmatrix
//endvimfold
	
//Uses the elements in the coulomn interaction matrix (HFdata), and the unitary matrix
//plus the orbital energies to construct a new unitary matric. Recoves small terms, and sorts
//the eigenenergies and corresponding eigenvectors, smallest elements first, for more stable 
//householder, QLDC (??) in the next iteration.
void HFmatrix::transform(double ** ddUnitaryMatrix){
//startvimfold

	//Temporary variables
	double Atemp;
	double Btemp;
	
	//First part; Coulomb energies. (Hartree and Fock energies)
	for (int alpha = 0; alpha < iNumberOfParticles; alpha++) {
		for (int gamma = 0; gamma < iNumberOfParticles; gamma++) {
			ppHFdata[alpha][gamma] = 0;
			for (int a = 0; a < iNumberOfParticles; a++) {
				for (int beta = 0; beta < iNumberOfParticles; beta++) {
					for (int delta = 0; delta < iNumberOfParticles; delta++) {
						//Need to ensure that 13 and 24 has same spin (correct acc. to the convention in this code).
						Atemp=0.0;
						Btemp=0.0;
						if ( alpha%2==delta%2 && beta%2==gamma%2 ){
						        Btemp = coulombIntegrals[alpha/2][beta/2][delta/2][gamma/2];
							}
						if ( alpha%2==gamma%2 && beta%2==delta%2 ){
								Atemp = coulombIntegrals[alpha/2][beta/2][gamma/2][delta/2];
						}						
						ppHFdata[alpha][gamma] += ddUnitaryMatrix[a][beta] * ddUnitaryMatrix[a][delta] * (Atemp-Btemp);
					}
            	}
			}
		}
	}
	
	//Second part: Adding eigenenergies of basis states.
	for (int i=0;i<iNumberOfParticles;i++){
		ppHFdata[i][i]+=pOrbitalEnergies[i/2];
	}

}//End of void HFmatrix::transform 
//endvimfold

//Changes rows of ddUnitaryMatrix to eigenvectors and pEigenvalues to the eigenvalues the HF-matrix
//Sorting eigenvalues in ddUnitaryMatrix
void HFmatrix::diagonalize(double* pEigenvalues, double** ddUnitaryMatrix){
//startvimfold

	//Temporary variable. Needed because of the way tred2 and tqli are programmed in lib.cpp
	double ddTrigonalMatrix[iNumberOfParticles];

	tred2(ppHFdata, iNumberOfParticles, pEigenvalues, ddTrigonalMatrix);
	tqli(pEigenvalues, ddTrigonalMatrix, iNumberOfParticles, ppHFdata);

	/*picksort:Num.Rec. inspired algo.
	O(n2). Not good for large iNumberOfParticles.
	Sorting dEigenvalues and pEigenvalues
	Smallest eigenvalues and corresponding eigenvectors first*/
	
	//Array containing indexes
	int iSortIndex[iNumberOfParticles];
	for (int i=0;i<iNumberOfParticles;i++){
		iSortIndex[i]=i;
	}

	int i,j;
	//Picksort algo
	double dEigTemp,dIndexTemp;
	for (i=iNumberOfParticles-2;i>=0;i--){
		dEigTemp=pEigenvalues[i];
		dIndexTemp=iSortIndex[i];//not neccesary
			for (j=i;j<=iNumberOfParticles-2;j++){
			pEigenvalues[j]=pEigenvalues[j+1];//cannot write i here
			iSortIndex[j]=iSortIndex[j+1];	
			if (pEigenvalues[j+1]>=dEigTemp) break;
		}
		pEigenvalues[j]=dEigTemp;
		iSortIndex[j]=dIndexTemp;//Can write i here
	}

	//Make unitary matrix from rows in orthonormal matrix
	for (i=0;i<iNumberOfParticles;i++){
		for (j=0;j<iNumberOfParticles;j++){
			ddUnitaryMatrix[i][j]=ppHFdata[i][iSortIndex[j]];
 		}
	}
	
	//Remove small elements that may lead to numerical errors
	for (i = 0; i < iNumberOfParticles; i++) {
		for (j = 0; j < iNumberOfParticles; j++) {
			if (fabs(ddUnitaryMatrix[i][j]) < 1e-4) {
				ddUnitaryMatrix[i][j] = 0.0;
			}
		}
	}
#if PRINT
	cout.precision(10);
	for (i=0;i<4;i++){
	cout<<pEigenvalues[i]<<","<<iSortIndex[i]<<"\n";
	}
cout<<"\n\n\n";
	for (i=0;i<4;i++){
	for (j=0;j<4;j++){
	cout<<ddUnitaryMatrix[i][j]<<"\t\t";
	}cout<<"\n";	}
cout<<"\n\n\n";
#endif

}//End of function HFmatrix::diagonalize: 
//endvimfold

//Finding <ab|V|cd>. This is the heaviest operation. Must be rewritten for larger other systems. 
//For large systems this operation should be done separately and stored. For basis functions with 3-dimensions
//MC-calculation of the integral may be the only practical solution.
double hFMatrixElements(int iNa,int iNb,int iNc,int iNd,int iNumberOfParticles, int iNrMeshpt, double dIntLimMin, double dIntLimMax){
//startvimfold

	double *W = new double[iNrMeshpt];
   	double *X = new double[iNrMeshpt];	

	int iZ = iNumberOfParticles; //Really electric charge of core
	
	//Kaller funksjon som returnerer meshpunkter X og vekter W
	gauleg(dIntLimMin,dIntLimMax,X,W,iNrMeshpt);
	
	//Summasjonsvar. Gir endelig svar
	double dInt_gl=0.0;
	
	//temporary variables
	double dInnerint,dInt1,dInt2;
	
	//Summation over all meshpoints
	for (int x_1=0;x_1<iNrMeshpt;x_1++){
		
		dInt1 = W[x_1] * X[x_1] *  radialWF(X[x_1],iNa,iZ) * radialWF(X[x_1],iNc,iZ);
		dInt2 = dInt1 * X[x_1];
		
		for (int x_2=0;x_2<iNrMeshpt;x_2++){
	
			//Computing the different factors in the integrand		
			//dInnerint = radialWF(X[x_2],iNb,iZ) * radialWF(X[x_2],iNd,iZ);
			
			//When x_1=x_2, the integrands are the same.		
			if (x_1>x_2) {
				//dInt_gl += dInt1 * W[x_2] * X[x_2] * X[x_2] * dInnerint;//
				dInt_gl += dInt1 * W[x_2] * X[x_2] * X[x_2] * radialWF(X[x_2],iNb,iZ) * radialWF(X[x_2],iNd,iZ);
			} else {
				//dInt_gl += dInt2 * W[x_2] * X[x_2] * dInnerint;//
				dInt_gl += dInt2 * W[x_2] * X[x_2] * radialWF(X[x_2],iNb,iZ) * radialWF(X[x_2],iNd,iZ);
			}
		}
	}
	
	//sletter arrays
	delete [] X;
	delete [] W;
	
	return dInt_gl;
}
//endvimfold

//Calculates radial part of the wavefunction.
double radialWF(double dR, int iN, int iZ){
//startvimfold

	double dZN =  (double)iZ/(double)iN;
	double dFac=1.0;

	//Calculate (N-1)!
	for (int  i = 2; i < iN; i++){
			dFac*=(double)i;
	}	
	return 	pow(2.0*dZN,1.5)*sqrt(dFac/(2*iN*iN*dFac))*
			exp(-dZN*dR)*
	 		LaguerreGeneral(iN - 1, 1, 2*dZN*dR );
}
//endvimfold

/*  Function to compute generalized Laguerre polynomials
    alpha cannot be smaller or equal than  -1.0
    The variable n has to be larger or equal 0
	*/
double LaguerreGeneral( int n, double alpha, double x){
//startvimfold
	double *glaguerre; 
	if ( alpha <= -1.0 ) {
	    cout << "LAGUERRE_GENERAL - Fatal error!" << endl;
	    cout << "The input value of ALPHA is=  "  << alpha  << endl;
	    cout << "but ALPHA must be greater than -1." << endl;
	}
	glaguerre = new double[n+1];
	if ( n >= 0 ) {
		for (int  i = 0; i <= n; i++) glaguerre[i] = 0.0;
			glaguerre[0] = 1.0;
		    if ( n > 0 ){
			      glaguerre[1] = 1.0+alpha-x;
		        // recursion relation for generalized Laguerre polynomials
		        for (int  i = 2; i <= n; i++){
	          		glaguerre[i] = ((2.0*i-1.0+alpha-x)*glaguerre[i-1]+
					(1.0-i-alpha)*glaguerre[i-2])/((float) i);
				}
			}
		}
	double GLaguerre = glaguerre[n];
	delete [] glaguerre;
	return GLaguerre;
}   // end function glaguerre
//endvimfold

int main(){
//startvimfold 

//Only needed in this method
	int iNumberOfParticles=4;

//declaring and constructing nxn matrix	
	double **ddUnitaryMatrix;
	ddUnitaryMatrix = (double **) matrix( iNumberOfParticles, iNumberOfParticles , sizeof(double));
/*for (int a=0; a< inumberofparticles; a++){
	ddunitarymatrix[a][a]=1.0;
	for (int b=a+1; b< inumberofparticles; b++){
 			ddunitarymatrix[a][b]=ddunitarymatrix[b][a]=0.0;
	}
}*/

	double* pEigenvalues = new double[iNumberOfParticles];
	double* pEigenvaluesOld = new double[iNumberOfParticles];// = new double[iNumberOfParticles];
			for (int i = 0; i < iNumberOfParticles; i++){
				pEigenvaluesOld[i] = 0; 	
			}

	bool bIterate = true;
	int iNumberOfIterations = 0;
	int iConverged;



	//Initializing hartree-fock matrix object
	HFmatrix A(2);
	
	//Start iterations
	while (bIterate) {
//		bIterate=false; //XXX XXX remove !!

//A.transform(ddUnity);
		A.diagonalize(pEigenvalues,ddUnitaryMatrix);
		A.transform(ddUnitaryMatrix);
//A.diagonalize(pEigenvalues,ddUnitaryMatrix);
//A.transform(ddUnitaryMatrix);
//A.diagonalize(pEigenvalues,ddUnitaryMatrix);
//A.transform(ddUnitaryMatrix);
//A.diagonalize(pEigenvalues,ddUnitaryMatrix);
#if 1
		if (iNumberOfIterations>0){
#endif
#if 1			
			iConverged = 0;
			for (int i = 0; i < iNumberOfParticles; i++) {
				if (abs((pEigenvaluesOld[i] - pEigenvalues[i])/pEigenvalues[i]) < 1e-6) { //init tol on top of program
					iConverged++;
				}
			}
			if (iConverged == iNumberOfParticles) {
				bIterate = false;
			}
			for (int i = 0; i < iNumberOfParticles; i++){
				pEigenvaluesOld[i] = pEigenvalues[i]; 	
			}
		}
#endif
		iNumberOfIterations++;
	}//End of while (bIterate)

cout << "Antall iterasjoner = " << iNumberOfIterations;

//Calculate energy
	//HFmatrix::<method> needed for calculating energy
	//method needed for writing data to file

//Write comments

//TESTING
/*for (int i=0;i<4;i++){
for (int j=0;j<4;j++){
cout<<", \t\t"<<ddUnitaryMatrix[i][j];
}
	cout<<"\n";
	}
*/
//delete[] pEigenvalues;
//delete[] pEigenvaluesOld;
}//End of Main().
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
