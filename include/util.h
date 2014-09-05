/**
 * \file
 * \brief Definition of the util class.
 */

#ifndef __UTIL_H__
#define __UTIL_H__ 1

//  includes from standard libraries
#include <iostream>
#include <vector>
#include <string>
//#include <ifstream>
#include <stdexcept>
#include <map>


namespace TCT {

  // \brief Util class for config file parsing

  class util {

    //private :

    public:

      std::map<std::string, std::string> _id_val;

      void parse(std::ifstream & cfgfile); 
  };
}
#endif
