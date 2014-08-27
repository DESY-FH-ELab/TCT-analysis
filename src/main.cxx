/**
 * \file
 * \brief Main program code.
 */

#include <iostream>
#include "sample.h"
#include "config.h"
#include "param.h"
#include "acquisition.h"
//#inculde "anaylyser.h" // inherits from class sample/measurement?

//using namespace TCT; // namespace of TCT_analysis is "TCT"
//using namespace std;

int main()
{
	std::cout << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;

	TCT::sample	dummyDUT;       // define DUT
	std::cout << dummyDUT << std::endl;	// print basic parameters of the sample


	std::string folder = "/home/hjansen/Sensors/testSensor";

	dummyDUT.ReadSampleCard(folder);	// read SampleCard and set parameters accordingly
	std::cout << dummyDUT << std::endl;

  	// analyser
 	//TCT::vdrift_ana vdrift;



	// print out method for analysers?!
	//std::cout << "dummy" << std::endl;

	std::cout << "end " << PACKAGE_NAME << std::endl;
  return 0;
}
