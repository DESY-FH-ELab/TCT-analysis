/**
 * \file
 * \brief Definition of the sample class.
 */

#ifndef __PULSE_H__
#define __PULSE_H__ 1

#include "laser.h"

// Stdandard includes
#include <iostream>
#include <vector>

namespace TCT {

	// \brief Abstract class for single pulse and avg pulse implementation (TCT and laser)

	class pulse_base {

	private :
		float _temp;	// in K
		float _volt;	// in Volt
		bool _IsEdgeTCT;
		double _x, _y; // do we need local and global coordinates ?? in um
		//string _sampleID; // not needed, sample will have a pulse variable 
		double SNR;
		//string _folder; 	// moved to sample
		//string _file;		// moved to sample
		double _AmpRefLaser;
	

	public :
		virtual double GetSNR() const = 0;
		void SetTemp( _temp = 0);
		
		laser type;
		laser settings;

	}; // end of class pulse

	// \brief Single pulse implementation

	class single_pulse : public pulse_base {

	public :
		virtual double GetSNR() const {
			return 1.0;
		}
	
	}; // end of Single pulse implementation

	class avg_pulse : public pulse_base {

	public :
		virtual double GetSNR() const {
			return 1.0;
		}
	
	}; // end of AVG pulse implementation

		
		

}
