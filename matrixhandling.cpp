#include "matrixhandling.h"

#include <complex>

#define d3 3

// Functions for d3xd3 matrices
void za(std::complex<double> c[d3][d3], std::complex<double> z, std::complex<double> a[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=z*a[i][j];
	}
}

void za(std::complex<double> c[d3][d3], double z, std::complex<double> a[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=std::complex<double>(z,0)*a[i][j];
	}
}

void aeb(std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		a[i][j]=b[i][j];
	}
}

void apb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=a[i][j]+c[i][j];
	}
}

void capb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=a[i][j]+b[i][j];
	}
}

void amb(std::complex<double> c[d3][d3],std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=a[i][j]-b[i][j];
	}
}
void axb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=0;
		for(int k=0;k<d3;k++)
			c[i][j]=c[i][j]+a[i][k]*b[k][j];
	}
}

void axbdag(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=0;
		for(int k=0;k<d3;k++)
			c[i][j]=c[i][j]+a[i][k]*conj(b[j][k]);
	}
}

void adagxb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=0;
		for(int k=0;k<d3;k++)
			c[i][j]=c[i][j]+conj(a[k][i])*b[k][j];
	}
}

void adagxbdag(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]){
	for(int i=0;i<d3;i++)
	for(int j=0;j<d3;j++){
		c[i][j]=0;
		for(int k=0;k<d3;k++)
			c[i][j]=c[i][j]+conj(a[k][i])*conj(b[j][k]);
	}
}

std::complex<double> multtrace(std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]){
	std::complex <double> tr;
	tr=std::complex<double>(0,0);

	for(int i=0;i<d3;i++)
	for(int k=0;k<d3;k++){
			tr=tr+a[i][k]*b[k][i];
	}
	return tr;
}

void adag(std::complex<double> a[d3][d3]){
       std::complex<double> tmp;
       for(int i=0;i<d3;i++){
               a[i][i]=conj(a[i][i]);
       }
       
       tmp=a[0][1];
       a[0][1]=conj(a[1][0]);
       a[1][0]=conj(tmp);

       tmp=a[0][2];
       a[0][2]=conj(a[2][0]);
       a[2][0]=conj(tmp);
       
       tmp=a[1][2];
       a[1][2]=conj(a[2][1]);
       a[2][1]=conj(tmp);
}

void expM(std::complex<double> U[d3][d3], std::complex<double> A[d3][d3]) {
  // U = exp(A)
  // U = 1 + sum_i 1/i! * A^i
  const int nmax=20;
  
  std::complex<double> Atmp1[d3][d3];
  std::complex<double> Atmp2[d3][d3];
  
  for(int i=0 ; i<d3; i++) {
    for(int j=0 ; j<d3; j++) {
      U[i][j] = 0;
    }
  }
  
  for(int i=0 ; i<d3; i++) {
    U[i][i] = 1;
    Atmp1[i][i] = 1;
  }

  double fact=1;

  for (int i=1; i<nmax; i++) {
    fact *= i;
    axb(Atmp2,Atmp1,A);
    aeb(Atmp1,Atmp2);
    za(Atmp2, 1.0/(double)fact, Atmp2);
    apb(U,Atmp2);
  }
}

void projA(std::complex<double> A[d3][d3], std::complex<double> U[d3][d3]) {
  // A = proj(U), where U is projected on a antihermitian traceless matrix A
  // proj(U) = ( U - U^dag)/2 - i*( U - U^dag)/2 )/3

  std::complex<double> Udag[d3][d3];
  std::complex<double> UmUdag[d3][d3];

  aeb(Udag,U);
  adag(Udag);

  amb(UmUdag, U,Udag);
  za(UmUdag, 0.5, UmUdag);

  aeb(A,UmUdag); 

  std::complex<double> trace;
  trace = UmUdag[0][0] + UmUdag[1][1] + UmUdag[2][2];
  for (int i=0; i<d3; i++) {
    A[i][i] -= trace/(double)3;
  }
}

double testAntiHerm(std::complex<double> A[d3][d3]) {
  // Define a norm for antihermitian properties
  // N_ah = || A || = | sum_i,j A_ij + conj(A_ji) |
  std::complex<double> cnorm;
  for (int i=0; i<d3; i++)
  for (int j=0; j<d3; j++) {
    cnorm += A[i][j] + conj(A[j][i]);
  }

  return abs(cnorm);
}

double testUnitarity(std::complex<double> U[d3][d3]) {
  // Test Unitarity U U^dag = 1
  // Define a norm for unitary properties
  // N_ah = || U || = ...

  std::complex<double> Udag[d3][d3];
  std::complex<double> UUdag[d3][d3];

  aeb(Udag,U);
  adag(Udag);
  axb(UUdag,U,Udag);

  std::complex<double> cnorm;
  for (int i=0; i<d3; i++) {
    for (int j=0; j<d3; j++) {
      cnorm += UUdag[i][j];
    }
    cnorm -= 1;
  }

  return abs(cnorm);
}
