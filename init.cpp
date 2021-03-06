/*
 * init.cc - program initialization
 *
 * Copyright © 2004 H.-P. Schadler  <hps@abyle.org>
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

#include "init.h"

#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>

#include "options.h"
#include "version.h"

int Init(int &argc, char *argv[], Options &opt) {
  const char texthelp[]="Usage: flow.x [OPTIONS] CONFIGFILE [TYPE]... \n"
      "Calculated the Wilson flow and writes measurements and gauge configurations\n"
      "as a function of Wilson flow time t.\n"
      "If -n/--nmeas > 1, the input file is interpreted as a list of input files!\n"
      "\n"
      "TYPE: b      Binary format (new) (default)\n"
      "      a      MILC format\n" 
      "      f      Fortran format\n"
      "\n"
      "Mandatory arguments to long options are mandatory for short options too.\n"
      "  -s, --Ns SSIZE         spatial lattice extent (default = 4)\n"
      "  -t, --Nt TSIZE         temporal lattice extent (default = 4)\n"
      "  -n, --nmeas NMEAS      number of configurations (default = 1)\n"
      "\n"  
      "  -e, --eps Epsilon      Flow time step\n"
      "  -f, --ftime t_max      Iterate Flow time until t=t_max\n"
      "\n"  
      "  -m  --meas             perform basic measurements and write them to file (default on)\n"
      "  -c  --writeconf        writes configurations to disk\n"
      "  -w  --wconfeveryt NUM   writes configurations to disk after every n-th t\n"
      "\n"  
      "  -d  --debug            write versbose informations to stdout\n"
      "\n"  
      "  -h  --help             display this help and exit\n"
      "  -v  --version          output version information and exit\n"
      "\n"
      "Exit status:\n"
      " 0  if OK,\n"
      " 1  if minor problems,\n"
      " 2  if serious trouble.\n"
      "\n"
      "Report bugs to hps@abyle.org";

  std::cout << "flow.x " << MAJOR_VERSION << "." << MINOR_VERSION << "." 
    << REVISION_VERSION << " ~ " << __DATE__ << " " << __TIME__ 
    << std::endl << std::endl
    << "Wilson flow calculations\n" 
    << "Calculated the Wilson flow and writes measurements and gauge configurations\n"
    << "as a function of Wilson flow time t.\n"
    << std::endl
    << "(C) 2014 Hans-Peter Schadler <hps@abyle.org>"
    << std::endl << std::endl;
  
  //if (argc<1) {
  //  std::cout << texthelp << std::endl;
  //  return 2;
  //}

  std::cout << std::endl << "Initializing... " << std::endl << std::endl;

  // Parse command line options and initialize global settings
  // We first set the default options for lattice size, measurements
  // and other important variables
  
  opt.ns=4;
  opt.nt=4;
  opt.nmeas=1;
  opt.meas=false;
  opt.writeconf=false;
  opt.eps=0.02;
  opt.tmax=2.00;
  opt.verbose=false;
  opt.iswriteceveryt=false;
  opt.writeceveryt=0;

  int c;

  while (1) {
    static struct option long_options[] =
    {
      /* These options don't set a flag.
      We distinguish them by their indices. */
      {"Ns", required_argument, 0, 's'},
      {"Nt", required_argument, 0, 't'},
      {"nmeas", required_argument, 0, 'n'},
      {"eps", required_argument, 0, 'e'},
      {"ftime", required_argument, 0, 'f'},
      /* These options set a flag. */
      {"meas", no_argument, 0, 'm'},
      {"writeconf", no_argument, 0, 'c'},
      {"wconfeveryt", required_argument, 0, 'w'},
      {"verbose", no_argument, 0, 'd'},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
		
	  /* getopt_long stores the option index here. */
	  int option_index = 0;

	  c = getopt_long (argc, argv, "s:t:n:e:f:w:mcdhv",
	  long_options, &option_index);

	  /* Detect the end of the options. */
	  if (c == -1) break;

	  switch (c) {
		  case 0:
			  /* If this option set a flag, do nothing else now. */
			  /*if (long_options[option_index].flag != 0)
				  break;
			  printf ("option %s", long_options[option_index].name);
			  if (optarg)
				  printf (" with arg %s", optarg);
			  printf ("\n");A */
	      //            if( strcmp( "3d", long_options[option_index].name ) == 0 )
		    //                 do3d = true;
			  break;

		  case 's':
		    opt.ns = std::atoi(optarg);
			  break;

		  case 't':
		    opt.nt = std::atoi(optarg);
			  break;
		  
      case 'n':
		    opt.nmeas = std::atoi(optarg);
			  break;
		  
      case 'e':
		    opt.eps = std::atof(optarg);
			  break;

		  case 'f':
		    opt.tmax = std::atof(optarg);
			  break;
      
      case 'w':
        opt.iswriteceveryt=true;
        opt.writeceveryt=std::atoi(optarg);
			  break;
			  
      case 'm':
        opt.meas = true;
        break;
      
      case 'c':
        opt.writeconf = true;
        break;
      
      case 'd':
        opt.verbose = true;
        break;

		  case 'v':
        std::cout << "flow.x " << MAJOR_VERSION << "." << MINOR_VERSION << "." 
          << REVISION_VERSION << " ~ " << __DATE__ << " " << __TIME__ 
          << std::endl << std::endl
          << "Wilson flow calculations" 
          << std::endl
          << "(C) 2014 Hans-Peter Schadler <hps@abyle.org>"
          << std::endl << std::endl;
			  exit(0);

		  case 'h':
			  std::cout << texthelp << std::endl;
			  exit(0);

		  default:
			  std::cout << texthelp << std::endl;
			  exit(0);
	  }
  }
  
  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    opt.filenamelist = argv[optind];
  }else{
    std::cout << "ERROR: No configuration file specified!" << std::endl;
    exit(1);
  }

  opt.filenames.resize(opt.nmeas);

  if (opt.nmeas>1) {
    // Read finname file and fill fevname
    std::ifstream fin;
    fin.open(opt.filenamelist.c_str());
    if (fin.is_open()!=true) {
      std::cout  << "ERROR: File " << opt.filenamelist <<  " to read configuration filename list could not be opened!" << std::endl;
      exit(1);
    }

    std::string strtmp;
    int n=0;
    while (n<opt.nmeas && getline(fin, opt.filenames.at(n))) {
      n++;
    }

    if (n<opt.nmeas) {
      std::cout << "WARNING: Only found " << n << " names in " << opt.filenamelist << " !" << std::endl;
      opt.nmeas=n;
    }
  } else {
    opt.filenames[0] = opt.filenamelist;
  }

  if (argc>optind+1) {
    if (argv[optind+1][0]=='b') {
      opt.type=0; // Binary2 (new format)
    }else if(argv[optind+1][0]=='a'){
      opt.type=1; // MILC
    }else if(argv[optind+1][0]=='f'){
      opt.type=2; // Fortran
    }else{
      std::cout << std::endl << "ERROR: Requested TYPE not understood!" << std::endl;
      std::cout << std::endl << texthelp << std::endl;
      return 1;
    }   
  }else{ // default is setting type to binary2
    opt.type=0;
  }

  return 0;
}
