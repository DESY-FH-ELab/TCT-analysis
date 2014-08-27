/**
 * \file
 * \brief Definition of the sample class.
 */

#ifndef __SAMPLE_H__
#define __SAMPLE_H__ 1

//  includes from standard library
#include <iostream>
//#include <vector>
#include<string>

//  includes from TCT classes
//#include "pulse.h"
//#include "laser.h"


namespace TCT {

	// \brief Sample class for sample implementation	
	
	class sample {

	private :
		std::string _SampleID;	// identifier for sample
		std::string _Folder;	// folder where data from the sample is stored
		//std::string _file;	// all measurement files belong to measurement class, not to sample class
		double _Thickness;
		void SetThickness(double thickness){ _Thickness = thickness;}


	public :

		// the empty constructor sets a dummy strong and a dummy thickness
		sample() : 
			_SampleID("DummyString"),
			_Thickness(-1)
		{}

		// Default copy constructer should be fine
		sample(const sample &)               = default;
		sample & operator = (const sample &) = default;

		// Dectructor
		~sample() = default;

		double & Thickness(){ return _Thickness;}
		const double & Thickness() const{ return _Thickness;}

		std::string & SampleID() {return _SampleID;}
		const std::string & SampleID() const {return _SampleID;}

		void ReadSampleCard(std::string folder);
	};

}

inline std::ostream & operator << (std::ostream & os, const TCT::sample & sam) {
	return os 	<< "\n This sample is called \n" << sam.SampleID()
			<< "   - thickness = " << sam.Thickness() << std::endl;
}
#endif
