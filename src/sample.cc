/**
 * \file
 * \brief Implementation of sample methods
 */

#include<string>

#include "sample.h"

namespace TCT {

  void sample::ReadSampleCard(){

    //read file from folder. fill name, thickness, effective doping etc

    std::cout << "\n reading SampleCard from " << _Folder << std::endl;
    double thick = 320.; // rad this from card
    SetThickness(thick);

    std::string sampleID = "DummySample1"; // read this from card
    SetSampleID(sampleID); 
  }


}
