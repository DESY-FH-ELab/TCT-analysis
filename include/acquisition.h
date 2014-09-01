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
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

namespace TCT {

  // \brief Abstract class for acquisition, acq_single and  acq_multi inherit from it

  class acquisition_base {

    private :
      float _Temp;	// in K
      //float _BiasVolt;	// in Volt
      bool _IsEdgeTCT;
      double _x, _y; // do we need local and global coordinates ?? in um
      double _SNR;
      double _AmpRefLaser;
      bool _IsManyPulseStructure;
      bool _IsMIP;
      uint32_t _Nsamples;
      float _SampleInterval;	// in ns

   protected:
      float _BiasVolt;	// in Volt

    public :
      acquisition_base() :
        _BiasVolt(500.)
      {};

      acquisition_base(float bias) :
        _BiasVolt(bias)
      {};

      virtual double GetSNR() const = 0;
      void SetTemp( float temp) { _Temp = 0;}

      uint32_t & Nsamples(){ return _Nsamples;}
      const uint32_t & Nsamples() const{ return _Nsamples;}

      float & SampleInterval(){ return _SampleInterval;}
      const float & SampleInterval() const{ return _SampleInterval;}

      std::vector<double> volt; // encapsulate ??
      std::vector<double> time;

      void SetBiasVolt(float volt) { _BiasVolt = volt;} // !! put back to base class
      float BiasVolt() {return _BiasVolt;} // !! put back to base class

  }; // end of class acquisition_base


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class acquisition_avg : public acquisition_base {

    private : 
      TProfile* _Profile;
      TProfile* _ProfileFILTERED;
      TH2F* _H_acqs2D;

    public :
      acquisition_avg() = default;//:
      //_IsMIP(false),
      //_Nsamples(1996),
      //_SampleInterval(0.1)
      //{}

      // Default copy constructer should be fine
      acquisition_avg(const acquisition_avg &)               = default;
      acquisition_avg & operator = (const acquisition_avg &) = default;

      // Dectructor
      ~acquisition_avg() = default;

      virtual double GetSNR() const {
	return 1.0;
      }
      TProfile* Profile() { return _Profile;}
      TProfile* ProfileFILTERED() { return _ProfileFILTERED;}
      TH2F* H_acqs() { return _H_acqs2D;}

  }; // end of acquisition_avg implementation


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // \brief Acquisition with single pulse implementation

  class acquisition_single : public acquisition_base {

    private :

      std::string _name;
      uint32_t ievt;
      float _maxamplitude;	// this is the amplitude of a single pulse, at MaxSigLocation
      float _amplneg;		// most neg value before rising edge
      float _amplpos;		// most pos value before rising edge
      float _amplneglate;	// most neg value after falling edge
      float _amplposlate;	// most neg value after falling edge
      float _avg;		// between edges for alpha 
      float _avgshort;		// avg shortly after risind edge (beta,p,...)
      float _delay;		// from acq start to start of RE, start of RE is defined by S2nRef in param.h
      float _delayfilt;		// from acq start to start of RE for filtered pulses
      float _width;		// from 50% of RE to 50% of FE
      float _rise;		// from start to position of maxamplitude 
      float _rise1090;		// rise time from 10% to 90%
      float _fall;		// from position of maxamp to end of pulse, end of FE is defined by S2nRef in param.h
      float _s2nval;		// avg / Noise
      float _Offset;		// baseline offset before RE
      float _Offset_end;	// baseline offset after FE
      float _Noise;		// Noise before RE
      float _Noise_end;	// Noise after FE
      TH1F *_H_acquisition;		// TH1F to store the data
      TH1F *_H_acquisitionFILTERED;	// TH1F to store the data, filtered
      uint32_t _Nsamples;
      float _SampleInterval;	// in ns


    public :
	acquisition_single() : 
	acquisition_base(500.),

      //acquisition_single() = default;
      acquisition_single(float bias) : 
	acquisition_base(bias),
      //  SetBiasVolt(500.);
      //};
	_SampleInterval(0.1){};
      //_Nsamples(1996),
      //_SampleInterval(0.1)


      // Default copy constructer should be fine
      acquisition_single(const acquisition_single &)               = default;
      acquisition_single & operator = (const acquisition_single &) = default;

      // Dectructor
      ~acquisition_single() = default;


      virtual double GetSNR() const {
	return 1.0;
      }

      //float BiasVolt() {return _BiasVolt;} // !! put back to base class

      void Print();
      void GetOffsetNoise(uint32_t iEv);
      bool Select();
      void DrawPulse();
      bool Read(FILE *infile, uint32_t iFile);
      bool Read(FILE *infile, uint32_t iFile, TCT::acquisition_avg avg); 
      void WriteNtuple();
      void SignalFinder();
      void SignalManipulator();
      void ClearStruct();
      void SetName(std::string name);
      void SetNsamples(uint32_t nsamples) { _Nsamples = nsamples;}
      uint32_t & Nsamples(){ return _Nsamples;}
      const uint32_t & Nsamples() const{ return _Nsamples;}

      float & SampleInterval(){ return _SampleInterval;}
      const float & SampleInterval() const{ return _SampleInterval;}

  }; // end of acquisition_single implementation

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class acquisition_multi_single : public acquisition_base {

    public :
      acquisition_multi_single() = default;//:
      //_IsMIP(false),
      //_Nsamples(1996),
      //_SampleInterval(0.1)
      //{}

      // Default copy constructer should be fine
      acquisition_multi_single(const acquisition_multi_single &)               = default;
      acquisition_multi_single & operator = (const acquisition_multi_single &) = default;

      // Dectructor
      ~acquisition_multi_single() = default;
      virtual double GetSNR() const {
	return 1.0;
      }

  }; // end of acquisition_multi_single implementation


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class acquisition_multi_avg : public acquisition_base {

    private : 
      TProfile* _Profile;

    public :
      acquisition_multi_avg() = default; //:
      //_IsMIP(false),
      //_Nsamples(1996),
      //_SampleInterval(0.1)
      //{}

      // Default copy constructer should be fine
      acquisition_multi_avg(const acquisition_multi_avg &)               = default;
      acquisition_multi_avg & operator = (const acquisition_multi_avg &) = default;

      // Dectructor
      ~acquisition_multi_avg() = default;

      virtual double GetSNR() const {
	return 1.0;
      }
      TProfile* Profile() {return _Profile;}

  }; // end of acquisition_multi_avg implementation	

}

#endif
