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
#include <fstream>
#include <stdexcept>
#include <map>


namespace TCT {

  // \brief Util class for config file parsing

  class util {

    private :

      std::map<std::string, std::string> _id_val;
      bool _IsRead = false;

    public:

      void parse(std::ifstream & cfgfile); 
      std::map<std::string, std::string> ID_val() {return _id_val;}

      bool IsRead() { return _IsRead;}

  };
}
#endif
