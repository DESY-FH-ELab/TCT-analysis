/**
 * \file
 * \brief Implementation of TCT::util methods
 */

//  includes from standard libraries
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <map>

//  includes from TCT classes
#include "util.h"

namespace TCT {

  void util::parse(std::ifstream & cfgfile) {
    std::string id, eq, val;
    _IsRead = false;

    while(cfgfile >> id ){

      if (id[0] == '#') {
        cfgfile.ignore(256,'\n');
        continue;  // skip comments
      }
      else if (id[0] == '[') {
        cfgfile.ignore(256,'\n');
        continue;  // skip group flag
      }
      else{ cfgfile >> eq >> val;

	if (eq != "=") throw std::runtime_error("Parse error");

	_id_val[id] = val;
      }

    }
#ifndef USE_GUI
    std::cout << "\n identifiers and values from passed file" << std::endl; 
    for(auto i : _id_val) {
      std::cout << i.first << " " << i.second << " " << "\n";
    }
#endif

    _IsRead = true;

    return;
  }
}
