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
      float _AmplPosEarly_Cut;	// in V
      float _AmplNegEarly_Cut;	// in V
      float _AmplPosLate_Cut;	// in V
      float _AmplNegLate_Cut;	// in V
      uint32_t _PrintEvent;

      std::string _OutFolder; // Folder, where the root file is writen to 
      std::string _OutSample_ID;
      std::string _OutTemp;
      std::string _OutVolt;
      std::string _OutPos; // for X and Y position



    public:

      analysis() :
	_Noise_Cut(0.01),
	_NoiseEnd_Cut(0.01),
        _S2n_Cut(9.),
	_S2n_Ref(2.),
	_Amplitude_Cut(0.01),
	_Width_Cut(3.0),
	_AmplPosEarly_Cut(0.02),
	_AmplNegEarly_Cut(-0.02),
	_AmplPosLate_Cut(0.01),
	_AmplNegLate_Cut(-0.02),
	_OutFolder("../results/"),
        _PrintEvent(-1)
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

      float AmplPosEarly_Cut() { return _AmplPosEarly_Cut;}
      void SetAmplPosEarly_Cut(float cut_value) { _AmplPosEarly_Cut = cut_value;}
      const float & AmplPosEarly_Cut() const { return _AmplPosEarly_Cut;}

      float AmplNegEarly_Cut() { return _AmplNegEarly_Cut;}
      void SetAmplNegEarly_Cut(float cut_value) { _AmplNegEarly_Cut = cut_value;}
      const float & AmplNegEarly_Cut() const { return _AmplNegEarly_Cut;}

      float AmplPosLate_Cut() { return _AmplPosLate_Cut;}
      void SetAmplPosLate_Cut(float cut_value) { _AmplPosLate_Cut = cut_value;}
      const float & AmplPosLate_Cut() const { return _AmplPosLate_Cut;}

      float AmplNegLate_Cut() { return _AmplNegLate_Cut;}
      void SetAmplNegLate_Cut(float cut_value) { _AmplNegLate_Cut = cut_value;}
      const float & AmplNegLate_Cut() const { return _AmplNegLate_Cut;}

      uint32_t PrintEvent() { return _PrintEvent;}
      void SetPrintEvent(int num ) { _PrintEvent = num;}
      const uint32_t & PrintEvent() const { return _PrintEvent;}

      //float () { return _;}
      //void Set(float ) { _ = ;}
      //const float & () const { return _;}

      void SetParameters(std::map<std::string, std::string> id_val);
      bool AcqsSelecter(TCT::acquisition_single *acq);
      void AcqsAnalyser(TCT::acquisition_single *acq, uint32_t iAcq, TCT::acquisition_avg *acqAvg);
      void AcqsProfileFiller(TCT::acquisition_single *acq, TCT::acquisition_avg *acqAvg);

      void AcqsWriter(std::string SampleID, std::string temp, std::string volt, std::vector<TCT::acquisition_single> *allAcqs, TCT::acquisition_avg *acqAvg);
      void AcqsWriter(std::string SampleID, std::string volt, std::vector<TCT::acquisition_single> *allAcqs, TCT::acquisition_avg *acqAvg);
      void AcqsWriter(std::vector<TCT::acquisition_single> *acqs, TCT::acquisition_avg *acqAvg);

      std::string OutFolder() {return _OutFolder;} 
      void SetOutFolder(std::string string) { _OutFolder = string;}
      const std::string & OutFolder() const {return _OutFolder;}
      
      std::string OutSample_ID() {return _OutSample_ID;} 
      void SetOutSample_ID(std::string string) { _OutSample_ID = string;}
      const std::string & OutSample_ID() const {return _OutSample_ID;}
      
      std::string OutTemp() {return _OutTemp;} 
      void SetOutTemp(std::string string) { _OutTemp = string;}
      const std::string & OutTemp() const {return _OutTemp;}
      
      std::string OutVolt() {return _OutVolt;} 
      void SetOutVolt(std::string string) { _OutVolt = string;}
      const std::string & OutVolt() const {return _OutVolt;}
      
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



