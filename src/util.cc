/**
 * \file
 * \brief Implementation of util methods
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

	//std::cout << "ID: " << id << "   eq: " << eq << "   val: " << val << std::endl;
	std::cout << "ID: " << id  << "   val: " << val << std::endl;
	if (eq != "=") throw std::runtime_error("Parse error");

	_id_val[id] = val;
      }

    }
    for(auto i : _id_val) {
      std::cout << i.first << " " << i.second << " " << "\n";
    }

    return;
  }
}
