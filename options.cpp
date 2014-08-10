/*
 * options.cc - global configuration
 *
 * Copyright © 2013 H.-P. Schadler  <hanspeter.schadler@uni-graz.at>
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

#include "options.h"

#include <iostream>
 
void Options::PrintSettings() {
  std::cout << "SETTINGS:" << std::endl
  << "Ns = " << ns << std::endl
  << "Nt = " << nt << std::endl << std::endl
  << "Nmeas = " << nmeas << std::endl;
  if (filenames.size() > 1) {
    std::cout << "Filename list = " << filenamelist << std::endl;
  }
  std::cout << "Filename = " << filenames[0] << std::endl;
  switch (type) {
    case 0:
      std::cout << "Storage format: binary (new)" << std::endl;
      break;

    case 1:
      std::cout << "Storage format: MILC text" << std::endl;
      break;

    case 2:
      std::cout << "Storage format: Fortran" << std::endl;
      break;

    default:
      std::cout << "Storage format: ERROR!!!" << std::endl;
  } 
  std:: cout << std::endl << "Write config. at t=t_max: " << writeconf << std::endl
  << "Write config. after every n-th t: " << writeceveryt << std::endl
  << "Write measurements: " << meas << std::endl << std::endl
  << "Wilson Flow settings:" << std::endl
  << "Flow step eps = " << eps << std::endl
  << "Max. flow time t_max = " << tmax << std::endl
  << "END SETTINGS" << std::endl;
}

