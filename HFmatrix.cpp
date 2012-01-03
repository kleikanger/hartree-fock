/*

balblabla

   */



#include <iostream>
#include "lib/lib.h"
#include "HFmatrix.h"
#include <cmath>

double calculateInnerIntegrals(double , int , int );  //REMOVE!

	
using std::cout;


//Not programmed by me.
double LaguerreGeneral(int,double,double);
//Programmed by me
double radialWF(double,int,int);
//double gaussianQuadrature(double,double,int,int,int,int,int,int); 
double hFMatrixElements(int ,int ,int ,int ,int ,int ,double, double);

//Constructs Hartre-Fock (nCutoff x nCutoffx) martix
HFmatrix::HFmatrix(int nCutoff){

	nCutoff=2;
	
	//XXX Should may be part of constructor call?
	iNumberOfParticles=2*nCutoff;
	iNrMeshpt=250;
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
	
	//Initialize elements of variable:coulombIntegrals
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
					cout<<"\n";
			    }
		 	}
		}

	cout<< "laguerre "<<LaguerreGeneral(2 - 1, 1, 2.0*4.0*  2.0  /((double) 2))<<"\t\t,";
	cout<< "radialwf" <<radialWF(2,1.,4)<<"\n";


	//Declaring and constructing unity matrix	
  	ddUnity = (double **) matrix( iNumberOfParticles, iNumberOfParticles , sizeof(double));
	for (int a=0; a< iNumberOfParticles; a++){
		ddUnity[a][a]=1.0;
		for (int b=a+1; b< iNumberOfParticles; b++){
 			ddUnity[a][b]=ddUnity[b][a]=0.0;
		}
	}
	
	//initializing eigenvalues(energies) of orbitals
	for (int i = 0; i < iNumberOfParticles; i++) {
		pOrbitalEnergies[i] = -((double) (iNumberOfParticles*iNumberOfParticles))/((double) (2*(i+1)*(i+1))); //iNumberOfParticles really Z!!
	}
    
	//Declare elements of variable:ppHFdata
 	ppHFdata= (double **) matrix( iNumberOfParticles, iNumberOfParticles, sizeof(double));

	//Finding initial matrix
	transform(ddUnity);
#if 0
	//Prin XXX test
	for (int i=0;i<4;i++){
		for (int j=0;j<4;j++){
		cout<<", \t"<<ppHFdata[i][j];
		}
	cout<<"\n";
	}

	for (int i=0;i<4;i++){
		for (int j=0;j<4;j++){
		cout<<", \t"<<ddUnity[i][j];
		}
	cout<<"\n";
	}
#endif
//XXX SLETT ARRAYS!

}//END of constructor HFmatrix
	
void HFmatrix::transform(double ** ddUnitaryMatrix){
	
	for (int alpha = 0; alpha < iNumberOfParticles; alpha++) {
		for (int gamma = 0; gamma < iNumberOfParticles; gamma++) {
    		// Bidraget fra Coulomb-matriseelementet:
			ppHFdata[alpha][gamma] = 0;
			for (int a = 0; a < iNumberOfParticles; a++) {
				for (int beta = 0; beta < iNumberOfParticles; beta++) {
					for (int delta = 0; delta < iNumberOfParticles; delta++) {
                		
						double Atemp=0.0;
						double Btemp=0.0;
						
						//XXX Following should be stated as a function later. will be used every iteration.
						//Need to check if states are perpendicular (because of different spinstates)
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
	//Adding eigenenergies
	for (int i=0;i<iNumberOfParticles;i++){
		ppHFdata[i][i]+=pOrbitalEnergies[i/2];
	}	
	


}//End of void HFmatrix::transform 



//Changes rows of ddUnitaryMatrix to eigenvectors and pEigenvalues to the eigenvalues the HF-matrix
//Sorting eigenvalues in ddUnitaryMatrix
void HFmatrix::diagonalize(double* pEigenvalues, double** ddUnitaryMatrix){

	//Should modify tred2,tqli instead of declareing this array
	double ddTrigonalMatrix[iNumberOfParticles];

	tred2(ppHFdata, iNumberOfParticles, pEigenvalues, ddTrigonalMatrix);
	tqli(pEigenvalues, ddTrigonalMatrix, iNumberOfParticles, ppHFdata);

	/*picksort:Num.Rec. inspired algo.
	Sorting dEigenvalues and pEigenvalues
	Smallest eigenvalues and corresponding eigenvectors first*/
	
	//Array with indexes
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
#if 1
	for (i=0;i<4;i++){
	cout<<pEigenvalues[i]<<","<<iSortIndex[i]<<"\n";
	}
cout<<"\n\n\n";
	for (i=0;i<4;i++){
	for (j=0;j<4;j++){
	cout<<ppHFdata[i][j]<<"\t\t";
	}cout<<"\n";	}
cout<<"\n\n\n";
#endif


}//End of function HFmatrix::diagonalize: 






/*Finds Hartree and Fock part of the  matrixelement iNa,iNb by integratian of the direct and exchange term.
  Remember that interchange of the terms a and b will give the same results. therefore only the elements of one trigonal
  part of the matrix will be neccasary to compute. This fact can be shown analytically (easily), but the calculation of the 
  full matrix may not give a perfect symmetric matrix(??), and then we will have problems with the diag. (Househoulder) routine.
 */


int main(){

int iNumberOfParticles=4;



//Declaring and constructing nxn matrix	
double **ddUnitaryMatrix;
ddUnitaryMatrix = (double **) matrix( iNumberOfParticles, iNumberOfParticles , sizeof(double));
/*for (int a=0; a< iNumberOfParticles; a++){
	ddUnitaryMatrix[a][a]=1.0;
	for (int b=a+1; b< iNumberOfParticles; b++){
 			ddUnitaryMatrix[a][b]=ddUnitaryMatrix[b][a]=0.0;
	}
}*/

double* pEigenvalues = new double[iNumberOfParticles];
double* pEigenvaluesOld = new double[iNumberOfParticles];



#if 0
if (iNumberOfIterations>1){
	for (int i = 0; i < iNumberOfParticles i++){
		for (int j = 0; j < iNumberOfParticles; j++){
			pEigenvaluesOld[i][j] = pEigenvalues[i][j] 	
		}
	}
}
#endif

HFmatrix A(2);
//A.transform(ddUnity);
A.diagonalize(pEigenvalues,ddUnitaryMatrix);

#if 0
if (iNumberOfIterations>1){
	for {int i = 0; i < iNumberOfParticles i++){
		for (int j = 0; j < iNumberOfParticles; j++){
			pEigenvaluesOld[i][j] = pEigenvalues[i][j] 	
		}
	}
}
#endif

#if 0
for (int nt nConverged = 0;
	for (int i = 0; i < 2*nCutoff; i++) {
		if (abs((dPrevious[i] - d[i])/d[i]) < tol) { //init tol on top of program
		nConverged++;
	}
}

if (nConverged == 2*nCutoff) {
	converged = true;
}

} //while loos stops here; while (convergence)
#endif

for (int i=0;i<4;i++){
for (int j=0;j<4;j++){
cout<<", \t\t"<<ddUnitaryMatrix[i][j];
}
	cout<<"\n";
	}
}//End of Main().


//Finding <ab|V|cd>
double hFMatrixElements(int iNa,int iNb,int iNc,int iNd,int iNumberOfParticles, int iNrMeshpt, double dIntLimMin, double dIntLimMax){
#if 1
	int n1,n2,n3,n4,Z,N; 
	n1=iNa;
	n2=iNb;
	n3=iNc;
	n4=iNd;
	Z=iNumberOfParticles;
	N=250.0;
	double intCutoff=250;

	    double *x = new double[N]; // Integrasjonspunkter
		double *w = new double[N]; // Integrasjonsvekter
			    
			    double integral = 0;
				    gauleg(0, intCutoff, x, w, N); // Finner x og w
					    for (int i = 0; i < N; i++) {
							        integral += w[i]*x[i]*x[i]*radialWF(x[i], n1,Z)*
										                    radialWF(x[i], n3,Z)*
															                    calculateInnerIntegrals(x[i], n2, n4);
									    }
						    
						    // Sletter:
						    delete[] x;
							    delete[] w;
								    
								    return integral;
}

/******************************************************************************/

// Funksjon som beregner det indre integralet i Coulomb-integralene ved bruk
// av Gauss-kvaderatur:
double calculateInnerIntegrals(double r1, int n1, int n2) {
int N=250;
double intCutoff=250;
int Z=4;
	    double *x = new double[N]; // Integrasjonspunkter
	    double *w = new double[N]; // Integrasjonsvekter
			    
	    // Første integral beregnes her (fra 0 til r1):
	    double integral1 = 0;
	    gauleg(0, r1, x, w, N); // Finner x og w
			    for (int i = 0; i < N; i++) {
				        integral1 += w[i]*x[i]*x[i]*radialWF(x[i], n1,Z)*
					                     radialWF(x[i], n2,Z);
			    }
			    integral1 /= r1;
							    
			    // Andre integral beregnes her (fra r1 til intCutoff):
			    double integral2 = 0;
			    gauleg(r1, intCutoff, x, w, N); // Finner x og w
				    for (int i = 0; i < N; i++) {
				        integral2 += w[i]*x[i]*radialWF(x[i], n1,Z)*
		                     radialWF(x[i], n2,Z);
													    }
										    
			    // Sletter.
			    delete[] x;
			    delete[] w;
								    
			    // Returnerer svaret:
			    return integral1 + integral2;
}
#else

	double *W = new double[iNrMeshpt];
   	double *X = new double[iNrMeshpt];	

	int iZ = iNumberOfParticles;
	
	//Kaller funksjon som returnerer meshpunkter X og vekter W
	gauleg(dIntLimMin,dIntLimMax,X,W,iNrMeshpt);
	
	//Summasjonsvar. Gir endelig svar
	double dInt_gl=0.0;
	
	//temporary variables
	double dInnerint,dInt1,dInt2;
	
	//Summation over all meshpoints
	for (int x_1=0;x_1<iNrMeshpt;x_1++){

		//Computing the different factors in the integrand		
		dInt1 = W[x_1] * X[x_1] *  radialWF(X[x_1],iNa,iZ) * radialWF(X[x_1],iNc,iZ);
		dInt2 = dInt1 * X[x_1];

		for (int x_2=0;x_2<iNrMeshpt;x_2++){
	
			//Computing the different factors in the integrand		
			dInnerint = radialWF(X[x_2],iNb,iZ) * radialWF(X[x_2],iNd,iZ);
			
			//When x_1=x_2, the integrands are the same.		
			if (X[x_1]<X[x_2]) {
				dInt_gl += dInt1 * W[x_2] * X[x_2] * X[x_2] * dInnerint;//
			} else {
				dInt_gl += dInt2 * W[x_2] * X[x_2] * dInnerint;//
			}
		}
	}
	
	//sletter arrays
	delete [] X;
	delete [] W;
	
	return dInt_gl;
} 
#endif

//Calculates radial part of the wavefunction.
double radialWF(double dR, int iN, int iZ){

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


/*  Function to compute generalized Laguerre polynomials
    alpha cannot be smaller or equal than  -1.0
    The variable n has to be larger or equal 0
	*/
double LaguerreGeneral( int n, double alpha, double x){
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



