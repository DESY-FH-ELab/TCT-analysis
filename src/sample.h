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
		std::string _sampleID;	// identifier for sample
		std::string _folder;	// folder where data from the sample is stored
		//std::string _file;	// all measurement files belong to measurement class, not to sample class
		double _thickness;
		void SetThickness(double thickness){ _thickness = thickness;}


	public :

		// the empty constructor sets a dummy strong and a dummy thickness
		sample() : 
		_sampleID("DummyString"),
		_thickness(-1){
		}

		// Default copy constructer should be fine
		sample(const sample &)               = default;
		sample & operator = (const sample &) = default;

		// Dectructor
		~sample() = default;

		double & GetThickness(){ return _thickness;}
		const double & GetThickness() const{ return _thickness;}

		std::string & sampleID() {return _sampleID;}
		const std::string & sampleID() const {return _sampleID;}

		void ReadSampleCard(std::string folder);
	};

}

inline std::ostream & operator << (std::ostream & os, const TCT::sample & sam) {
	return os 	<< "\n This sample is called " << sam.sampleID()
			<< "\n   - thickness = " << sam.GetThickness() << std::endl;
}
#endif
