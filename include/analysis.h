/**
 * \file
 * \brief Definition of the analysis class.
 */

#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__ 1

// Standard includes
#include <iostream>
#include <vector>
#include <map>

// TCT includes
#include "acquisition.h"
#include "sample.h"

// ROOT includes

namespace TCT {

  // \brief analysis class for setting analysis-related parameters

  class analysis {

    private :

      float _Noise_Cut;		// in V
      float _NoiseEnd_Cut;
      float _S2n_Cut;
      float _S2n_Ref;
      float _Amplitude_Cut;	// in V
      float _Width_Cut;	// in ns

      std::string _OutFolder; // Folder, where the root file is writen to !! move to Analysis class



    public:

      analysis() :
	_Noise_Cut(0.01),
	_NoiseEnd_Cut(0.01),
        _S2n_Cut(9.),
	_S2n_Ref(2.),
	_Amplitude_Cut(0.01),
	_Width_Cut(3.0),
	_OutFolder("../results/")
      {
        //std::cout << "\n   *** No parameter map passes, using default cut values! ***" << std::endl;
      }


      analysis(std::map<std::string, std::string> id_val) :
        analysis()
      {
        SetParameters(id_val);
      }


      float S2n_Cut() { return _S2n_Cut;}
      void SetS2n_Cut (float val) { _S2n_Cut = val;}
      const float & S2n_Cut() const { return _S2n_Cut;}

      float S2n_Ref() { return _S2n_Ref;}
      void SetS2n_Ref (float val) { _S2n_Ref = val;}
      const float & S2n_Ref() const { return _S2n_Ref;}

      float Noise_Cut() { return _Noise_Cut;}
      void SetNoise_Cut(float cut_value) { _Noise_Cut = cut_value;}

      float NoiseEnd_Cut() { return _NoiseEnd_Cut;}
      void SetNoiseEnd_Cut(float cut_value) { _NoiseEnd_Cut = cut_value;}

      float Width_Cut() { return _Width_Cut;}
      void SetWidth_Cut(float cut_value) { _Width_Cut = cut_value;}
      const float & Width_Cut() const { return _Width_Cut;}

      float Amplitude_Cut() { return _Amplitude_Cut;}
      void SetAmplitude_Cut(float cut_value) { _Amplitude_Cut = cut_value;}
      const float & Amplitude_Cut() const { return _Amplitude_Cut;}

      void SetParameters(std::map<std::string, std::string> id_val);
      bool AcqsSelecter(TCT::acquisition_single *acq);
      void AcqsAnalyser(TCT::acquisition_single *acq, uint32_t iAcq, TCT::acquisition_avg *acqAvg);
      void AcqsProfileFiller(TCT::acquisition_single *acq, TCT::acquisition_avg *acqAvg);
      void AcqsWriter(TCT::sample* sample, std::vector<TCT::acquisition_single> *acqs, TCT::acquisition_avg *acqAvg);

      std::string OutFolder() {return _OutFolder;} // !! move to ana class
      void SetOutFolder(std::string string) { _OutFolder = string;} // !! move to ana class
      const std::string & OutFolder() const {return _OutFolder;}
      
    };

}

inline std::ostream & operator << (std::ostream & os, const TCT::analysis & ana) {
  return os	<< "\n   Parameters of this analysis are" 
  		<< "\n    S2n_Cut = " << ana.S2n_Cut() 
  		<< "\n    S2n_Ref = " << ana.S2n_Ref() 
  		<< "\n    Width_Cut = " << ana.Width_Cut() 
    << std::endl;
}

#endif



