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

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "ConfigData.hpp"
#include "init.h"
#include "matrixhandling.h"
#include "measuretime.h"
#include "options.h"

Options opt;

void PrintMatrix(std::complex<double> M[3][3]) {
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

std::complex<double> CalcPoll(ConfigData *config) {
        // Calculates the Polyakov loop spatial average
        std::complex<double> poll, trace;
        poll=std::complex<double> (0,0);
        trace=std::complex<double> (0,0);

        std::complex<double> up[3][3], uu[3][3], upaux[3][3];

        int is0=0;
        int i4=0;

        for (int i1=0; i1<opt.ns; i1++) {
                for (int i2=0; i2<opt.ns; i2++) {
                        for (int i3=0; i3<opt.ns; i3++) {
                                i4 = 0;
                                is0 = i1 + i2*opt.ns + i3*opt.ns*opt.ns + i4*opt.ns*opt.ns*opt.ns;

                                config->extract(*up, 3, is0);

                                for (i4=1; i4<opt.nt-1; i4++) {
                                        is0 = i1 + i2*opt.ns + i3*opt.ns*opt.ns + i4*opt.ns*opt.ns*opt.ns;
                                        config->extract(*uu, 3, is0);
                                        axb(upaux, up, uu);
                                        aeb(up, upaux);
                                }

                                i4 = opt.nt-1;
                                is0 = i1 + i2*opt.ns + i3*opt.ns*opt.ns + i4*opt.ns*opt.ns*opt.ns;

                                config->extract(*uu, 3, is0);
                                trace=multtrace(up, uu);

                                poll += trace;
                        }
                }
        }

        poll = poll/((double)3*opt.ns*opt.ns*opt.ns);

        return poll;
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

std::complex<double> CalcE(ConfigData *config) {
  // Calculates the plaquette spatial average over all 6N plaquettes
  int ispmu, ispnu;
  std::complex<double> sumplaqs = std::complex<double>(0,0), trace;
  std::complex<double> u1[3][3], u2[3][3], u3[3][3], u4[3][3], u23[3][3], u234[3][3];

  for (int is=0; is<opt.ns*opt.ns*opt.ns*opt.nt; is++) {
          for (int imu=0; imu<4; imu++) {
                  for (int inu=imu+1; inu<4; inu++) {
                          ispmu = config->neib(is, imu);
                          ispnu = config->neib(is, inu);

                          config->extract(*u1, imu, is);
                          config->extract(*u2, inu, ispmu);
                          config->extract(*u3, imu, ispnu);
                          config->extract(*u4, inu, is);

                          axbdag(u23, u2, u3);
                          axbdag(u234, u23, u4);
                          trace = multtrace(u1, u234);

                          sumplaqs += trace;
                  }
          }
  }

  return 1.0 - sumplaqs/(double)(6*3*opt.ns*opt.ns*opt.ns*opt.nt);
}

std::complex<double> CalcPlaq(ConfigData *config) {
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

  // Additional factor of 4 because of double (quad) counting of plaquettes
  return plaq/(double)(4*6*3*opt.ns*opt.ns*opt.ns*opt.nt);
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
  
void SmallFlowStep(ConfigData *config) {
  // Here we have to implement the RK scheme from 1006.4518 [hep-lat]
  // Appendix C
  // W_0 = V_t
  // W_1 = exp(1/4 Z_0) W_0
  // W_2 = exp(8/9 Z_1 - 17/26 Z_0) W_1
  // W_t+eps = exp(3/4 Z_2 - 8/9 Z_1 + 17/36 Z_0) W_2
  //
  // Z_i = eps Z(W_i)
  //
  // The input configuration 'config' will be overwritten with the
  // configuration at flow time t'=t+eps.
  
  std::complex<double> U[3][3];
  std::complex<double> S[3][3];
  std::complex<double> US[3][3];
  std::complex<double> A[3][3];
  std::complex<double> A2[3][3];
  std::complex<double> expA[3][3];
  std::complex<double> newU[3][3];

  // We allocate a new ConfigData
  // W_0 = V_t : W0 = config;
  ConfigData *W0 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  CopyConfig(W0, config);
  
  // W_1 = exp(1/4 Z_0) W_0
  // Multiplicate every link U_x,mu in W0 with 
  // exp(1/4 eps * U_x,mu*S_x,mu) where S_x,mu is the sum of staples around U_x,mu
  // and save it in W1.
  ConfigData *W1 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
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
  
  // W_2 = exp(8/9 Z_1 - 17/26 Z_0) W_1
  // Multiplicate every link U_x,mu in W_1 (config1) with 
  // exp(8/9 eps * V_x,mu*T_x,mu - 17/36 * U_x,mu*S_x,mu ) where T_x,mu \in W_1 and
  // S_x,mu \in W_0 and save it in W_2 (config2).
  ConfigData *W2 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
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

  // W_2 = exp(8/9 Z_1 - 17/26 Z_0) W_1
  // Multiplicate every link U_x,mu in W_1 (config1) with 
  // exp(8/9 eps * V_x,mu*T_x,mu - 17/36 * U_x,mu*S_x,mu ) where T_x,mu \in W_1 and
  // S_x,mu \in W_0 and save it in W_2 (config2).
  ConfigData *Wtpeps = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
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

  // Make Wtpeps the new configuration at t' = t+eps
  CopyConfig(config, Wtpeps);
  
  // Cleaning up
  delete W0; W0 = 0;
  delete W1; W1 = 0;
  delete W2; W2 = 0;
  delete Wtpeps; Wtpeps = 0;
}

int main(int argc, char *argv[]) {
  double tstart, tend;

  // Init program and parse commandline parameters
  Init(argc, argv, opt);
  opt.PrintSettings();

  std::cout << std::endl;

  ConfigData *config;

  // Read the config given on command line
  config = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  int ret=config->readBinaryConfig2(opt.filename);
  // int ret=config->MILCreadConfig(opt.filename);
  if (ret!=0) {
    std::cout << "ERROR: Reading MILC config file!" << std::endl;
    return 1;
  }
  
  // Test of staple sum
  std::complex<double> plaq;
  plaq = CalcPlaq(config);
  std::cout << "Plaquette (from staples) = " <<  plaq << std::endl;

  std::ofstream file;
  file.open("plaq.data");
  file << "t <Plaq_t> <E_t> <Poll_t>" << std::endl;
  
  std::cout << std::endl << "Flow calculation: ( t <Plaq_t> <E_t> <Poll_t> )" << std::endl;

  std::complex<double> E, poll;

  // Step t=0
  E = CalcE(config);
  poll = CalcPoll(config);
  std::cout << 0 << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;
  file << 0 << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;

  // Evolution in t
  for(double t=0; t<4; t+=opt.eps) {
    tstart = gettime();
    SmallFlowStep(config);
    plaq = CalcPlaq(config);
    E = CalcE(config);
    poll = CalcPoll(config);
    std::cout << t+opt.eps << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;
    file << t+opt.eps << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;
    tend = gettime();
    std::cout << "Evolution step done in " << tend - tstart << " s." << std::endl;

  }

  file.close();
}
