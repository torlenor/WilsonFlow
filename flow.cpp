/*
 * flow.cpp - Wilson flow calculations
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

#include <string>

#include "ConfigData.hpp"
#include "init.h"
#include "matrixhandling.h"
#include "options.h"

Options opt;

ConfigData *config;

void PrintMatrix(std::complex<double> M[3][3]) {
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void CalcStapleSum(std::complex<double> S[3][3], int x, int mu) {
  std::complex<double> U1[3][3];
  std::complex<double> U2[3][3];
  std::complex<double> U3[3][3];
  std::complex<double> U21[3][3];
  std::complex<double> U321[3][3];

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) {
    S[i][j]=std::complex<double>(0,0);
  }
  
  int xpmu = config->neib(x,mu); // Always the same

  for (int nu=0; nu<4; nu++) {
    if (mu != nu) {
      // S_x,mu = U_x+mu,nu * U^dag_x+nu,mu * U^dag_x,nu + U^dag_x-nu+mu,nu * U^dag_x-nu,mu * U_x-nu,nu
    
      // Upper part of staple
      int xpnu = config->neib(x,nu);

      config->extract(*U1, nu, x); // U_x,nu
      config->extract(*U2, mu, xpnu); // U_x+nu,mu
      config->extract(*U3, nu, xpmu);

      adagxbdag(U21, U2, U1);
      axb(U321, U3, U21);
      apb(S, U321);
      
      // Lower part of staple
      int xmnu = config->neib(x,nu+4);
      int xmnupmu = config->neib(xmnu,mu);

      config->extract(*U1, nu, xmnu); // U_x-nu,nu
      config->extract(*U2, mu, xmnu); // U_x-nu,mu
      config->extract(*U3, nu, xmnupmu); // U_x-nu+mu,nu

      adagxb(U21, U2, U1);
      adagxb(U321, U3, U21);
      apb(S, U321);
    }
  }
}

std::complex<double> CalcPlaq() {
  std::complex<double> U[3][3];
  std::complex<double> S[3][3];
  
  std::complex<double> plaq;
  for (int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for (int mu=0; mu<4; mu++) {
      config->extract(*U, mu ,x);
      CalcStapleSum(S, x, mu);
      plaq += multtrace(U,S);
    }
  }
  
  return plaq/(double)(6*3*opt.ns*opt.ns*opt.ns*opt.nt)/(double)4;
}

int main(int argc, char *argv[]) {
  // Init program and parse commandline parameters
  Init(argc, argv, opt);
  opt.PrintSettings();

  std::cout << std::endl;

  // Read the config given on command line
  config = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  int ret=config->MILCreadConfig(opt.filename);
  if (ret!=0) {
    std::cout << "ERROR: Reading MILC config file!" << std::endl;
    return 1;
  }
  
  // Test of staple sum
  std::complex<double> plaq;
  plaq = CalcPlaq();
  std::cout << "Plaquette (from staples) = " <<  plaq << std::endl;

  std::complex<double> I(0,1);

  std::complex<double> S[3][3] = {I,2,I,I,1,2,1,I,2};
  int x=10, mu=0;
  CalcStapleSum(S, x, mu);

  PrintMatrix(S);
  std::cout << testAntiHerm(S) << std::endl;
  std::cout << std::endl;

  projA(S,S);

  PrintMatrix(S);
  std::cout << testAntiHerm(S) << std::endl;
  std::cout << std::endl;
  
  std::complex<double> U[3][3];
  expM(U,S);
  
  PrintMatrix(U);
  std::cout << testUnitarity(U) << std::endl;
  std::cout << std::endl;

}
