/**
 * \file
 * \brief Implementation of util methods
 */

//  includes from standard libraries
#include <iostream>
#include <vector>
#include <string>
#include <ifstream>
#include <stdexcept>
#include <map>

//  includes from TCT classes
#include "util.h"


namespace TCT {

  void parse(std::ifstream & cfgfile) {
    std::string id, eq, val;

    while(cfgfile >> id >> eq >> val)
    {
      if (id[0] == '#') continue;  // skip comments
      if (eq != "=") throw std::runtime_error("Parse error");

      _id_val[id] = val;
    }
  }
