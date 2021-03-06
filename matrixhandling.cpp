/*
 * matrixhandling.cpp - Matrix functions
 * Last changes: 2014-08-09 - Cleanups
 *
 * Copyright © 2014 H.-P. Schadler  <hps@abyle.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "matrixhandling.h"

#include <complex>

#define d3 3

// Functions for d3xd3 matrices
void za(std::complex<double> c[d3][d3], std::complex<double> z, std::complex<double> a[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = z*a[i][j];
	}
}

void za(std::complex<double> c[d3][d3], double z, std::complex<double> a[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = std::complex<double>(z,0)*a[i][j];
	}
}

void aeb(std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]){
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		a[i][j] = b[i][j];
	}
}

void apb (std::complex<double> c[d3][d3], std::complex<double> a[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = a[i][j]+c[i][j];
	}
}

void capb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = a[i][j]+b[i][j];
	}
}

void amb(std::complex<double> c[d3][d3],std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = a[i][j]-b[i][j];
	}
}
void axb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = 0;
		for (int k=0; k<d3; k++)
			c[i][j] = c[i][j]+a[i][k]*b[k][j];
	}
}

void axbdag(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = 0;
		for (int k=0; k<d3; k++)
			c[i][j] = c[i][j]+a[i][k]*conj(b[j][k]);
	}
}

void adagxb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = 0;
		for (int k=0; k<d3; k++)
			c[i][j] = c[i][j]+conj(a[k][i])*b[k][j];
	}
}

void adagxbdag(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3],std::complex<double> b[d3][d3]) {
	for (int i=0; i<d3; i++)
	for (int j=0; j<d3; j++) {
		c[i][j] = 0;
		for (int k=0; k<d3; k++)
			c[i][j] = c[i][j]+conj(a[k][i])*conj(b[j][k]);
	}
}

std::complex<double> multtrace(std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]) {
	std::complex <double> tr(0,0);

	for (int i=0; i<d3; i++)
	for (int k=0; k<d3; k++) {
			tr = tr+a[i][k]*b[k][i];
	}
	return tr;
}

void adag(std::complex<double> a[d3][d3]) {
       std::complex<double> tmp;
       for (int i=0; i<d3; i++) {
        a[i][i] = conj(a[i][i]);
       }
       
       tmp = a[0][1];
       a[0][1] = conj(a[1][0]);
       a[1][0] = conj(tmp);

       tmp = a[0][2];
       a[0][2] = conj(a[2][0]);
       a[2][0] = conj(tmp);
       
       tmp = a[1][2];
       a[1][2] = conj(a[2][1]);
       a[2][1] = conj(tmp);
}

inline bool IsSmall(std::complex<double> U[d3][d3]) {
  double abs=0;
  for (int j=0; j<d3; j++)
  for (int i=0; i<d3; i++) {
    abs += std::abs(U[j][i]);
  }
  if (abs < 1e-15) {
    return true;
  } else {
    return false;
  }
}

void expM(std::complex<double> U[d3][d3], std::complex<double> A[d3][d3]) {
  // U = exp(A)
  // U = 1 + sum_i 1/i! * A^i
  const int nmax=20;
  
  std::complex<double> Atmp1[d3][d3], Atmp2[d3][d3];
  
  for (int i=0; i<d3; i++) {
    for (int j=0; j<d3; j++) {
      U[i][j] = 0;
    }
  }
  
  for (int i=0; i<d3; i++) {
    U[i][i] = 1;
    Atmp1[i][i] = 1;
  }

  double fact=1;

  int i=0;
  for (i=1; i<nmax; i++) {
    fact *= i;
    axb(Atmp2,Atmp1,A);
    aeb(Atmp1,Atmp2);
    za(Atmp2, 1.0/(double)fact, Atmp2);
    apb(U,Atmp2);
    if (IsSmall(Atmp2))
      break;
  }
}

void projA(std::complex<double> A[d3][d3], std::complex<double> U[d3][d3]) {
  // A = proj(U), where U is projected on a antihermitian traceless matrix A
  // proj(U) = ( U - U^dag)/2 - i*( U - U^dag)/2 )/3

  std::complex<double> Udag[d3][d3], UmUdag[d3][d3];

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

  std::complex<double> Udag[d3][d3], UUdag[d3][d3];

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
