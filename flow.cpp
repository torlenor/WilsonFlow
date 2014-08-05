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
// #include "jackknife.h"
#include "matrixhandling.h"
#include "measuretime.h"
#include "options.h"
#include "su3.h"

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

        config->extract(*up, is0, 3);

        for (i4=1; i4<opt.nt-1; i4++) {
          is0 = i1 + i2*opt.ns + i3*opt.ns*opt.ns + i4*opt.ns*opt.ns*opt.ns;
          config->extract(*uu, is0, 3);
          axb(upaux, up, uu);
          aeb(up, upaux);
        }

        i4 = opt.nt-1;
        is0 = i1 + i2*opt.ns + i3*opt.ns*opt.ns + i4*opt.ns*opt.ns*opt.ns;

        config->extract(*uu, is0, 3);
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

      config->extract(*U1, x, nu); // U_x,nu
      config->extract(*U2, xpnu, mu); // U_x+nu,mu
      config->extract(*U3, xpmu, nu);

      adagxbdag(U21, U2, U1);
      axb(U321, U3, U21);
      apb(S, U321);
      
      // Lower part of staple
      int xmnu = config->neib(x,nu+4);
      int xmnupmu = config->neib(xmnu,mu);

      config->extract(*U1, xmnu, nu); // U_x-nu,nu
      config->extract(*U2, xmnu, mu); // U_x-nu,mu
      config->extract(*U3, xmnupmu, nu); // U_x-nu+mu,nu

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

        config->extract(*u1, is, imu);
        config->extract(*u2, ispmu, inu);
        config->extract(*u3, ispnu, imu);
        config->extract(*u4, is, inu);

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
      config->extract(*U, x, mu);
      CalcStapleSum(S, config, x, mu);
      plaq += multtrace(U,S);
    }
  }

  // Additional factor of 4 because of double (quad) counting of plaquettes
  return plaq/(double)(4*6*3*opt.ns*opt.ns*opt.ns*opt.nt);
}

template <class T>
void Jackknife(std::vector<T> &in, T &mean, T &error) {
  int nmeas=in.size();

  mean = 0;
  for (int j=0; j<nmeas; j++) {
    mean += in[j];
  }   
  mean=mean/(double)nmeas;

  error = 0;
  for (int j=0; j<nmeas; j++){
    error += pow(mean-in[j],2);
  }   

  error = sqrt((double)(nmeas-1)*error/(double)(nmeas));
}

template <class T>
void CalcJackMean(std::vector<T> &meas, T &mean, T &error) {
  std::vector<T> jackdata;
  jackdata.resize(meas.size());

  for (int n=0; n<meas.size(); n++) {
    jackdata[n] = 0;
    for (int j=0; j<meas.size(); j++) {
      if (n!=j)
        jackdata[n] += meas[j];
    }   
    jackdata[n] = jackdata[n]/(double)(opt.nmeas-1);
  }

  Jackknife(jackdata, mean, error);
}

void CopyConfig(ConfigData *configDst, ConfigData *configSrc) {
  std::complex<double> U[3][3];

  for(int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for(int mu=0; mu<4; mu++) {
      configSrc->extract(*U, x, mu);
      configDst->replace(*U, x, mu);
    }
  }
}

void CalcZ(std::complex<double> Z[3][3], ConfigData *W, int x, int mu) {
  std::complex<double> U[3][3];
  std::complex<double> S[3][3];
  std::complex<double> US[3][3];

  W->extract(*U, x, mu);
  CalcStapleSum(S, W, x, mu);
  axb(US, U, S);
  projA(Z,US);
  za(Z,-1.0,Z);
}
  
void SmallFlowStep(ConfigData *config) {
  // Here we have to implement the RK scheme from 1006.4518 [hep-lat]
  // Appendix C
  // W_0 = V_t
  // W_1 = exp(1/4 Z_0) W_0
  // W_2 = exp(8/9 Z_1 - 17/36 Z_0) W_1
  // W_t+eps = exp(3/4 Z_2 - 8/9 Z_1 + 17/36 Z_0) W_2
  //
  // Z_i = eps Z(W_i)
  //
  // The input configuration 'config' will be overwritten with the
  // configuration at flow time t'=t+eps.
  
  std::complex<double> U[3][3];

  std::complex<double> Z0[3][3];
  std::complex<double> Z1[3][3];
  std::complex<double> Z2[3][3];
  
  std::complex<double> A[3][3];
  std::complex<double> A2[3][3];

  std::complex<double> expA[3][3];
  std::complex<double> newU[3][3];

  // We allocate a new ConfigData
  // W_0 = V_t : W0 = config;
  ConfigData *W0 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  ConfigData *S0 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  CopyConfig(W0, config);
  
  // W_1 = exp(1/4 Z_0) W_0
  ConfigData *W1 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  for (int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for (int mu=0; mu<4; mu++) {
      // For link U_x,mu
      CalcZ(Z0, W0, x, mu);
      S0->replace(*Z0, x, mu); // save this into ConfigData for Z0s
      za(A, 1.0/4.0*opt.eps, Z0);

      expM(expA, A);
      W0->extract(*U, x, mu);
      axb(newU, expA, U);

      if (testUnitarity(newU) > 1e-5) {
        std::cout << "WARNING: New link in SmallFlowStep() not unitary anymore!" << std::endl;
      }

      W1->replace(*newU, x, mu);
    }
  }
  
  // W_2 = exp(8/9 Z_1 - 17/36 Z_0) W_1
  ConfigData *W2 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  ConfigData *S1 = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  for (int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for (int mu=0; mu<4; mu++) {
      // First part with Z0
      // CalcZ(Z0, W0, x, mu);
      S0->extract(*Z0, x, mu); // save this into ConfigData for Z0s
      za(A, -17.0/36.0*opt.eps, Z0);

      // Add second part with Z1
      CalcZ(Z1, W1, x, mu);
      S1->replace(*Z1, x, mu); // save this into ConfigData for Z0s
      za(A2, 8.0/9.0*opt.eps, Z1);

      // Add both parts
      apb(A,A2);

      expM(expA, A);
      W1->extract(*U, x, mu);
      axb(newU, expA, U);

      if (testUnitarity(newU) > 1e-5) {
        std::cout << "WARNING: New link in SmallFlowStep() not unitary anymore!" << std::endl;
      }

      W2->replace(*newU, x, mu);
    }
  }

  // W_t+eps = exp(3/4 Z_2 - 8/9 Z_1 + 17/36 Z_0) W_2
  ConfigData *Wtpeps = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
  for (int x=0; x<opt.ns*opt.ns*opt.ns*opt.nt; x++) {
    for (int mu=0; mu<4; mu++) {
      // For link U_x,mu
      
      // First part with Z0
      S0->extract(*Z0, x, mu); // load Z0
      za(A, 17.0/36.0*opt.eps, Z0);

      // Add second part with Z1
      S1->extract(*Z1, x, mu); // load Z1
      za(A2, -8.0/9.0*opt.eps, Z1);
      
      // Add both parts
      apb(A,A2);
      
      // Add third part with Z2
      CalcZ(Z2, W2, x, mu);
      za(A2, 3.0/4.0*opt.eps, Z2);
      
      // Add both parts
      apb(A,A2);

      expM(expA, A);
      W2->extract(*U, x, mu);
      axb(newU, expA, U);

      if (testUnitarity(newU) > 1e-5) {
        std::cout << "WARNING: New link in SmallFlowStep() not unitary anymore!" << std::endl;
      }

      Wtpeps->replace(*newU, x, mu);
    }
  }

  // Make Wtpeps the current configuration now at t' = t+eps
  CopyConfig(config, Wtpeps);
  
  // Cleaning up
  delete W0; W0 = 0;
  delete W1; W1 = 0;
  delete W2; W2 = 0;
  delete S0; S0 = 0;
  delete S1; S1 = 0;
  delete Wtpeps; Wtpeps = 0;
}

int main(int argc, char *argv[]) {
  double tstart, tend;

  // Init program and parse commandline parameters
  Init(argc, argv, opt);
  opt.PrintSettings();

  std::cout << std::endl;

  std::vector<double> plaqdata, edata, polldata;

  // Start meas loop over gauge configuraions
  ConfigData *config;
  for (int n=0; n<opt.nmeas; n++) {
    std::cout << std::endl << "Measurement " << n+1 << " of " << opt.nmeas << std::endl;

    // Read gauge configuration from file
    config = new ConfigData(opt.ns, opt.ns, opt.ns, opt.nt, 3);
    int ret;
    switch (opt.type) {
      case 0:
        ret=config->readBinaryConfig2(opt.filenames[n]);
        if(ret!=0){
          std::cout << "ERROR: Reading binary config file (new storage format)!" << std::endl;
          return 1;
        }   
        break;

      case 1:
        ret=config->MILCreadConfig(opt.filenames[n]);
        if(ret!=0){
          std::cout << "ERROR: Reading MILC config file!" << std::endl;
          return 1;
        }   
        break;

      case 2:
        ret=config->readFConfig(opt.filenames[n]);
        if(ret!=0){
          std::cout << "ERROR: Reading Fortran config file!" << std::endl;
          return 1;
        }
        break;

      default:
        std::cout << "ERROR: Requested TYPE not understood!" << std::endl;
        return 1;
      }   
    
    // Test of staple sum
    std::complex<double> plaq;
    plaq = CalcPlaq(config);
    std::cout << "Plaquette (from staples) = " <<  plaq << std::endl;

    std::stringstream fmeasname;
    fmeasname << opt.filenames[n] << ".flow";
    std::ofstream file;
    file.open(fmeasname.str().c_str());
    if (! file.is_open() ) {
      std::cout << "There was a problem opening the file " << fmeasname.str() << " !" << std::endl;
      exit(1);
    }
    file << "# Input: " << opt.filenames[n] << std::endl;
    file << "# t <Plaq_t> <E_t> <Poll_t>" << std::endl;
    
    if (opt.verbose)
      std::cout << std::endl << "Flow calculation: ( t <Plaq_t> <E_t> <Poll_t> )" << std::endl;

    std::complex<double> E, poll;

    // Step t=0
    // E = CalcE(config);
    E = 1.0 - plaq;
    poll = CalcPoll(config);

    if (opt.verbose)
      std::cout << 0 << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;

    file << 0 << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;

    // Evolution in t up to t_max in steps of eps
    for (double t=0; t<opt.tmax; t+=opt.eps) {
      std::cout << "Flow time t = " << t+opt.eps << " ... " << std::flush;
      tstart = gettime();
      SmallFlowStep(config);
      plaq = CalcPlaq(config);
      // E = CalcE(config);
      E = 1.0 - plaq;
      poll = CalcPoll(config);

      if (opt.verbose)
        std::cout << t+opt.eps << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;

      file << t+opt.eps << " " << std::real(plaq) << " " << std::real(E) << " " << std::abs(poll) << std::endl;
      tend = gettime();
      std::cout << " done in " << tend - tstart << " s." << std::endl;
    }

    file.close();
    
    plaqdata.push_back(std::real(plaq));
    polldata.push_back(std::abs(poll));
    edata.push_back(std::real(E));

    delete config; config=0; // cleanup
  }
    
  // Calculate mean values using Jackknife
  // void CalcJackMean(std::vector<T> &meas, T &mean, T &error)
  double plaqmean, plaqerr, pollmean, pollerr, emean, eerr;
  CalcJackMean(plaqdata, plaqmean, plaqerr);
  CalcJackMean(polldata, pollmean, pollerr);
  CalcJackMean(edata, emean, eerr);

  // Write mean values calculated using Jackknife
  std::stringstream fmeanname;
  fmeanname << "mean_" << opt.ns << "x" << opt.nt << ".flow";
  std::ofstream mfile;
  mfile.open(fmeanname.str().c_str());
  
  if (! mfile.is_open() ) {
    std::cout << "There was a problem opening the file " << fmeanname.str() << " !" << std::endl;
    exit(1);
  }

  mfile << "# t = " << opt.tmax << std::endl;
  mfile << "# <Plaq_t> <Plaq_t>err <E_t> <E_t>err <Poll_t> <Poll_t>err" << std::endl;
  mfile << plaqmean << " " << plaqerr << " "
        << emean << " " << eerr << " "
        << pollmean << " " << pollerr << " "
        << std::endl;

  mfile.close();
}
