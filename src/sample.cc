/**
 * \file
 * \brief Implementation of sample methods
 */

#include<string>

#include "sample.h"

namespace TCT {

	void sample::ReadSampleCard(std::string folder){
	
	//read file from folder, fill thickness, effective doping etc
		
		std::cout << "reading SampleCard from " << folder << std::endl;
		_Thickness = 500.0;
	}
}
