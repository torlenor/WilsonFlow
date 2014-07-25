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

#include <iostream>
#include <fstream>
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

void CalcStapleSum(std::complex<double> S[3][3], ConfigData *config, int x, int mu) {
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

std::complex<double> CalcE() {
  // Calculates the plaquette spatial average over all 6N plaquettes
  int ispmu, ispnu;
  std::complex<double> sumplaqs = std::complex<double>(0,0), trace;
  std::complex<double> u1[3][3], u2[3][3],u3[3][3], u4[3][3], u23[3][3], u234[3][3];

  for(int is = 0;is<opt.ns*opt.ns*opt.ns*opt.nt;is++){
          for(int imu = 0;imu<4;imu++){
                  for(int inu = imu+1;inu<4;inu++){
                          ispmu = config->neib(imu, is);
                          ispnu = config->neib(inu, is);

                          config->extract(*u1,imu,is);
                          config->extract(*u2,inu,ispmu);
                          config->extract(*u3,imu,ispnu);
                          config->extract(*u4,inu,is);

                          axbdag(u23,u2,u3);
                          axbdag(u234,u23,u4);
                          trace=multtrace(u1,u234);

                          sumplaqs = sumplaqs + trace;
                  }
          }
  }

  sumplaqs=sumplaqs/((double)6*3*(double)opt.ns*opt.ns*opt.ns*opt.nt); // factor because sum over N lattice points and md from trace
  return sumplaqs;

  // return (double)3*opt.ns*opt.ns*opt.ns*opt.nt - sumplaqs;
}

std::complex<double> CalcPlaq() {
  std::complex<double> U[3][3];
  std::complex<double> S[3][3];
  
  std::complex<double> plaq;
  for (int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for (int mu=0; mu<4; mu++) {
      config->extract(*U, mu ,x);
      CalcStapleSum(S, config, x, mu);
      plaq += multtrace(U,S);
    }
  }
  
  return plaq/(double)(6*3*opt.ns*opt.ns*opt.ns*opt.nt)/(double)4;
}

void CopyConfig(ConfigData *configDst, ConfigData *configSrc) {
  std::complex<double> U[3][3];

  for(int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for(int mu=0; mu<4; mu++) {
      configSrc->extract(*U, mu, x);
      configDst->replace(*U, mu, x);
    }
  }

}

void SmallFlowStep() {
  // Here we have to implement the RK scheme from 1006.4518 [hep-lat]
  // Appendix C
  // W_0 = V_t
  // W_1 = exp(1/4 Z_0) W_0
  // W_2 = exp(8/9 Z_1 - 17/26 Z_0) W_1
  // W_t+eps = exp(3/4 Z_2 - 8/9 Z_1 + 17/36 Z_0) W_2
  //
  // Z_i = eps Z(W_i)
  
  std::complex<double> U[3][3];
  std::complex<double> S[3][3];
  std::complex<double> US[3][3];
  std::complex<double> A[3][3];
  std::complex<double> A2[3][3];
  std::complex<double> expA[3][3];
  std::complex<double> newU[3][3];

  // We allocate a new ConfigData
  ConfigData *W0 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  // W_0 = V_t : W0 = config;
  CopyConfig(W0, config);
  
  ConfigData *W1 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  // W_1 = exp(1/4 Z_0) W_0
  // Multiplicate every link U_x,mu in W0 with 
  // exp(1/4 eps * U_x,mu*S_x,mu) where S_x,mu is the sum of staples around U_x,mu
  // and save it in W1.
  for(int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for(int mu=0; mu<4; mu++) {
      // For link U_x,mu
      W0->extract(*U, mu, x);
      CalcStapleSum(S, W0, x, mu);
      axb(US, U, S);
      projA(A,US);
      za(A,1.0/4.0*opt.eps,A);
      expM(expA, A);
      axb(newU, expA, U);
      if (testUnitarity(newU) > 1e-5) {
        std::cout << "WARNING: New link in SmallFlowStep() not unitary anymore!" << std::endl;
      }
      W1->replace(*newU, mu, x);
    }
  }
  
  ConfigData *W2 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  // W_2 = exp(8/9 Z_1 - 17/26 Z_0) W_1
  // Multiplicate every link U_x,mu in W_1 (config1) with 
  // exp(8/9 eps * V_x,mu*T_x,mu - 17/36 * U_x,mu*S_x,mu ) where T_x,mu \in W_1 and
  // S_x,mu \in W_0 and save it in W_2 (config2).
  for(int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for(int mu=0; mu<4; mu++) {
      // For link U_x,mu
      
      // First part with Z0
      W0->extract(*U, mu, x);
      CalcStapleSum(S, W0, x, mu);
      axb(US, U, S);
      projA(A,US);
      za(A,-17.0/36.0*opt.eps,A);

      // Add second part with Z1
      W1->extract(*U, mu, x);
      CalcStapleSum(S, W1, x, mu);
      axb(US, U, S);
      projA(A2,US);
      za(A2,8.0/9.0*opt.eps,A2);

      // Add both parts
      apb(A,A2);

      expM(expA, A);
      axb(newU, expA, U);
      if (testUnitarity(newU) > 1e-5) {
        std::cout << "WARNING: New link in SmallFlowStep() not unitary anymore!" << std::endl;
      }
      W2->replace(*newU, mu, x);
    }
  }

  ConfigData *Wtpeps = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  // W_2 = exp(8/9 Z_1 - 17/26 Z_0) W_1
  // Multiplicate every link U_x,mu in W_1 (config1) with 
  // exp(8/9 eps * V_x,mu*T_x,mu - 17/36 * U_x,mu*S_x,mu ) where T_x,mu \in W_1 and
  // S_x,mu \in W_0 and save it in W_2 (config2).
  for(int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for(int mu=0; mu<4; mu++) {
      // For link U_x,mu
      
      // First part with Z0
      W0->extract(*U, mu, x);
      CalcStapleSum(S, W0, x, mu);
      axb(US, U, S);
      projA(A,US);
      za(A,17.0/36.0*opt.eps,A);

      // Add second part with Z1
      W1->extract(*U, mu, x);
      CalcStapleSum(S, W1, x, mu);
      axb(US, U, S);
      projA(A2,US);
      za(A2,-8.0/9.0*opt.eps,A2);
      
      // Add both parts
      apb(A,A2);
      
      // Add third part with Z2
      W2->extract(*U, mu, x);
      CalcStapleSum(S, W1, x, mu);
      axb(US, U, S);
      projA(A2,US);
      za(A2,3.0/4.0*opt.eps,A2);
      
      // Add both parts
      apb(A,A2);

      expM(expA, A);
      axb(newU, expA, U);
      if (testUnitarity(newU) > 1e-5) {
        std::cout << "WARNING: New link in SmallFlowStep() not unitary anymore!" << std::endl;
      }
      Wtpeps->replace(*newU, mu, x);
    }
  }

  CopyConfig(config, Wtpeps);
  
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
  CalcStapleSum(S, config, x, mu);

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

  std::ofstream file;
  file.open("plaq.data");

  std::complex<double> E;
  for(double t=0; t<4; t+=opt.eps) {
    SmallFlowStep();
    plaq = CalcPlaq();
    E = CalcE();
    std::cout << t << " " << plaq << std::endl;
    file << t << " " << std::real(plaq) << " " << std::real(E) << std::endl;
  }
}
