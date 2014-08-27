/**
 * \file
 * \brief Definition of the sample class.
 */

#ifndef __ACQUISITION_H__
#define __ACQUISITION_H__ 1

//#include "laser.h" // not needed, acq doesnt need to know frmo which laser it was produced

// Standard includes
#include <iostream>
#include <vector>
#include "/home/hjansen/root/include/TH1F.h"

namespace TCT {

	// \brief Abstract class for acquisition, acq_single and  acq_multi inherit from it

	class pulse_base {

	private :
		float _Temp;	// in K
		float _BiasVolt;	// in Volt
		bool _IsEdgeTCT;
		double _x, _y; // do we need local and global coordinates ?? in um
		//string _sampleID; // not needed, sample will have a pulse variable 
		double _SNR;
		double _AmpRefLaser;
		bool IsManyPulseStructure;
	

	public :
		virtual double GetSNR() const = 0;
		void SetTemp( float temp) { _Temp = 0;}
		
	}; // end of class pulse

	
	// \brief Single pulse implementation

	class pulse_single : public pulse_base {

	private :

		std::vector<double> volt;
		std::vector<double> time;
		char* Name[20];
		uint32_t ievt;
		float maxamplitude;	// this is the amplitude of a single pulse, at MaxSigLocation
		float amplneg;		// most neg value before rising edge
		float amplpos;		// most pos value before rising edge
		float amplneglate;	// most neg value after falling edge
		float amplposlate;	// most neg value after falling edge
		float avg;		// between edges for alpha 
		float avgshort;		// avg shortly after risind edge (beta,p,...)
		float delay;		// from acq start to start of RE, start of RE is defined by S2nRef in param.h
		float delayfilt;	// from acq start to start of RE for filtered pulses
		float width;		// from 50% of RE to 50% of FE
		float rise;		// from start to position of maxamplitude 
		float rise1090;		// rise time from 10% to 90%
		float fall;		// from position of maxamp to end of pulse, end of FE is defined by S2nRef in param.h
		float s2nval;		// avg / Noise
		float Offset;		// baseline offset before RE
		float Offset_end;	// baseline offset after FE
		float Noise;		// Noise before RE
		float Noise_end;	// Noise after FE
		TH1F *Hpulse;		// TH1F to store the data
		TH1F *HpulseFILTERED;	// TH1F to store the data, filtered


	public :
		virtual double GetSNR() const {
			return 1.0;
		}


		void Print();
		void GetOffsetNoise(uint32_t iEv);
		bool Select();
		void DrawPulse();
		bool Read(FILE *infile, uint32_t iFile);
		void WriteNtuple();
		void SignalFinder();
		void SignalManipulator();
		void ClearStruct();
	
	}; // end of Single pulse implementation


	class pulse_avg : public pulse_base {

	public :
		virtual double GetSNR() const {
			return 1.0;
		}
	
	}; // end of pulse_avg implementation


	class pulse_multi : public pulse_base {

	public :
		virtual double GetSNR() const {
			return 1.0;
		}
	
	}; // end of pulse_multi implementation

		
		

}

#endif
