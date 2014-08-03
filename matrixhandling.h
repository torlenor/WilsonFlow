#ifndef _MATRIXHANDLING_H
#define _MATRIXHANDLING_H

#include <complex>

#define d3 3

// Functions for d3xd3 matrices
void aeb(std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]);
void apb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3]);
void capb(std::complex<double> c[d3][d3], std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]);
void amb(std::complex<double> c[d3][d3],std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]);

void axb(std::complex<double> a[d3][d3],std::complex<double> b[d3][d3], std::complex<double> c[d3][d3]);
void axbdag(std::complex<double> a[d3][d3],std::complex<double> b[d3][d3], std::complex<double> c[d3][d3]);
void adagxb(std::complex<double> a[d3][d3],std::complex<double> b[d3][d3], std::complex<double> c[d3][d3]);
void adagxbdag(std::complex<double> a[d3][d3],std::complex<double> b[d3][d3], std::complex<double> c[d3][d3]);
std::complex<double> multtrace(std::complex<double> a[d3][d3], std::complex<double> b[d3][d3]);
void za(std::complex<double> c[d3][d3], double z, std::complex<double> a[d3][d3]);
void adag(std::complex<double> a[d3][d3]);

void expM(std::complex<double> U[d3][d3], std::complex<double> A[d3][d3]);
void projA(std::complex<double> A[d3][d3], std::complex<double> U[d3][d3]);

double testAntiHerm(std::complex<double> A[d3][d3]);
double testUnitarity(std::complex<double> U[d3][d3]);

#endif /* _MATRIXHANDLING_H */
