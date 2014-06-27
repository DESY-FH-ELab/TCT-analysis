/**
 * \file
 * \brief Definition of the sample class.
 */

#ifndef __SAMPLE_H__
#define __SAMPLE_H__ 1

#include "pulse.h"
#include "laser.h"

// Stdandard includes
#include <iostream>
#include <vector>

namespace TCT {

	// \brief Sample class for sample implementation	
	
	class sample {

	private :
		double _thickness;
		string _sampleID;
		string _folder;
		string _file;

	};

}
