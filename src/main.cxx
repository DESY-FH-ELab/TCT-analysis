/**
 * \file
 * \brief Main program code.
 */

#include <iostream>
#include "sample.h"
#include "config.h"
//#inculde "anaylyser.h" // class or methods of sample ?

//using namespace TCT; // namespace of TCT_analysis is "TCT"
//using namespace std;

int main()
{
	TCT::sample	dummyDUT;       // define DUT
	std::cout << dummyDUT << std::endl;	// print basic parameters of the sample

	std::cout << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;

	std::string folder = "/home/hjansen/Sensors/testSensor";

	dummyDUT.ReadSampleCard(folder);	// read SampleCard and set parameters accordingly

  	// analyser
 	//TCT::vdrift_ana vdrift;



	// print out method for analysers?!
	//std::cout << "dummy" << std::endl;

  return 0;
}
