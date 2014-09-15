/**
 * \file
 * \brief Implementation of sample methods
 */

//  STD includes
#include<string>

//  TCT includes
#include "sample.h"

namespace TCT {

  void sample::SetParameters(std::map<std::string, std::string> id_val){

#ifdef DEBUG 
    std::cout << " start SAMPLE::SetParameterd: Read and set parameters from input map" << std::endl; 
#endif

    for( auto i : id_val){
      //if(i.first == "") _ = atof((i.second).c_str());
      if(i.first == "SampleCardFolder") _Folder = i.second;
      if(i.first == "SampleID") _SampleID = i.second;
      if(i.first == "Thickness") _Thickness = atof((i.second).c_str());
    }

#ifdef DEBUG 
    std::cout << " end SAMPLE::SetParameters" << std::endl; 
#endif

    return;
  }

}
