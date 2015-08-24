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

      uint32_t _MaxAcqs;
      uint32_t _Mode;
      uint32_t _TCT_Mode;
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

      float _SampleThickness;

      std::string _OutFolder; // Folder, where the root file is writen to 
      std::string _OutSample_ID;
      std::string _OutSubFolder;
      std::string _OutSubsubFolder;
      std::string _OutPos; // for X and Y position

      bool _DoSmearing;
      float _AddNoise;
      float _AddJitter;

      std::string _SampleCard; // 
      std::string _DataFolder; // 
      bool _SaveToFile;
      bool _SaveSingles;
      bool _LeCroyRAW;

      //scanning parameters
      uint32_t _CH_Det;
      uint32_t _CH_PhDiode;
      uint32_t _CH_Trig;
      uint32_t _OptAxis;
      bool _DO_focus;
      uint32_t _FPerp;
      float _FFWHM;
      float _FTlow;
      float _FThigh;
      float _FDLow;
      float _FDHigh;
      bool _FSeparateCharges;
      float _Movements_dt;

    public:

      analysis() :
	_Mode(0),
    _MaxAcqs(-1),
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
        _PrintEvent(-1),
	_DoSmearing(false),
	_AddNoise(.0),
	_AddJitter(.0),
	_SampleCard("def"),
	_DataFolder("def"),
	_SaveToFile(false),
	_SaveSingles(false),
    _LeCroyRAW(false),
    _CH_Det(0),
    _CH_PhDiode(0),
    _CH_Trig(0),
    _OptAxis(3),
    _DO_focus(false),
    _Movements_dt(0),
    _TCT_Mode(0)
      {
        //std::cout << "\n   *** No parameter map passes, using default cut values! ***" << std::endl;
      }


      analysis(std::map<std::string, std::string> id_val) :
        analysis()
      {
        SetParameters(id_val);
      }

      uint32_t Mode() { return _Mode;}
      void SetMode(uint32_t val) { _Mode = val;}
      const uint32_t & Mode() const { return _Mode;}

      uint32_t TCT_Mode() { return _TCT_Mode;}
      void SetTCT_Mode(uint32_t val) { _TCT_Mode = val;}
      const uint32_t & TCT_Mode() const { return _TCT_Mode;}

      uint32_t MaxAcqs() { return _MaxAcqs;}
      void SetMaxAcqs(uint32_t val) { _MaxAcqs = val;}
      const uint32_t & MaxAcqs() const { return _MaxAcqs;}

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

      bool DoSmearing() { return _DoSmearing;}
      void SetDoSmearing(bool val) { _DoSmearing = val;}
      const bool & DoSmearing() const { return _DoSmearing;}

      bool SaveSingles() { return _SaveSingles;}
      void SetSaveSingles(bool val) { _SaveSingles = val;}
      const bool & SaveSingles() const { return _SaveSingles;}

      bool LeCroyRAW() { return _LeCroyRAW;}
      void SetLeCroyRAW(bool val) { _LeCroyRAW = val;}
      const bool & LeCroyRAW() const { return _LeCroyRAW;}

      bool SaveToFile() { return _SaveToFile;}
      void SetSaveToFile(bool val) { _SaveToFile = val;}
      const bool & SaveToFile() const { return _SaveToFile;}

      float AddNoise() { return _AddNoise;}
      void SetAddNoise(float noise) { _AddNoise = noise;}
      const float & AddNoise() const { return _AddNoise;}

      float AddJitter() { return _AddJitter;}
      void SetAddJitter(float val) { _AddJitter = val;}
      const float & AddJitter() const { return _AddJitter;}

      //begin scanning section
      uint32_t CH_Det() { return _CH_Det;}
      void SetCH_Det(uint32_t val) { _CH_Det = val;}
      const uint32_t & CH_Det() const { return _CH_Det;}

      uint32_t CH_PhDiode() { return _CH_PhDiode;}
      void SetCH_PhDiode(uint32_t val) { _CH_PhDiode = val;}
      const uint32_t & CH_PhDiode() const { return _CH_PhDiode;}

      uint32_t CH_Trig() { return _CH_Trig;}
      void SetCH_Trig(uint32_t val) { _CH_Trig = val;}
      const uint32_t & CH_Trig() const { return _CH_Trig;}

      uint32_t OptAxis() { return _OptAxis;}
      void SetOptAxis(uint32_t val) { _OptAxis = val;}
      const uint32_t & OptAxis() const { return _OptAxis;}

      bool DO_focus() { return _DO_focus;}
      void SetDO_focus(bool val) { _DO_focus = val;}
      const bool & DO_focus() const { return _DO_focus;}

      uint32_t FPerp() { return _FPerp;}
      void SetFPerp(uint32_t val) { _FPerp = val;}
      const uint32_t & FPerp() const { return _FPerp;}

      float FFWHM() { return _FFWHM;}
      void SetFFWHM(float val) { _FFWHM = val;}
      const float & FFWHM() const { return _FFWHM;}

      float FTlow() { return _FTlow;}
      void SetFTlow(float val) { _FTlow = val;}
      const float & FTlow() const { return _FTlow;}

      float FThigh() { return _FThigh;}
      void SetFThigh(float val) { _FThigh = val;}
      const float & FThigh() const { return _FThigh;}

      float FDLow() { return _FDLow;}
      void SetFDLow(float val) { _FDLow = val;}
      const float & FDLow() const { return _FDLow;}

      float FDHigh() { return _FDHigh;}
      void SetFDHigh(float val) { _FDHigh = val;}
      const float & FDHigh() const { return _FDHigh;}

      bool FSeparateCharges() { return _FSeparateCharges;}
      void SetFSeparateCharges(bool val) { _FSeparateCharges = val;}
      const bool & FSeparateCharges() const { return _FSeparateCharges;}

      float Movements_dt() { return _Movements_dt;}
      void SetMovements_dt(float val) { _Movements_dt = val;}
      const float & Movements_dt() const { return _Movements_dt;}

      //end scanning section

      std::string SampleCard () { return _SampleCard;}
      void SetSampleCard(std::string val) { _SampleCard = val;}
      const std::string & SampleCard() const { return _SampleCard;}

      std::string DataFolder () { return _DataFolder;}
      void SetDataFolder(std::string val) { _DataFolder = val;}
      const std::string & DataFolder() const { return _DataFolder;}

      //float () { return _;}
      //void Set(float ) { _ = ;}
      //const float & () const { return _;}

      void SetParameters(std::map<std::string, std::string> id_val);
      bool AcqsSelecter(TCT::acquisition_single *acq);
      void AcqsSmearer(TCT::acquisition_single *acq, float noise, bool);
      void AcqsSmearer(TCT::acquisition_single *acq, bool, float jitter);
      void AcqsAnalyser(TCT::acquisition_single *acq, uint32_t iAcq, TCT::acquisition_avg *acqAvg);
      void AcqsProfileFiller(TCT::acquisition_single *acq, TCT::acquisition_avg *acqAvg);

      void AcqsWriter(std::vector<TCT::acquisition_single> *acqs, TCT::acquisition_avg *acqAvg, bool HasSubs = true);
      //void AcqsWriterNoSubs(std::vector<TCT::acquisition_single> *acqs, TCT::acquisition_avg *acqAvg);

      std::string OutFolder() {return _OutFolder;} 
      void SetOutFolder(std::string string) { _OutFolder = string;}
      const std::string & OutFolder() const {return _OutFolder;}
      
      std::string OutSample_ID() {return _OutSample_ID;} 
      void SetOutSample_ID(std::string string) { _OutSample_ID = string;}
      const std::string & OutSample_ID() const {return _OutSample_ID;}

      float SampleThickness() { return _SampleThickness;}
      void SetSampleThickness(float val) { _SampleThickness = val;}
      const float & SampleThickness() const { return _SampleThickness;}
      
      std::string OutSubFolder() {return _OutSubFolder;} 
      void SetOutSubFolder(std::string string) { _OutSubFolder = string;}
      const std::string & OutSubFolder() const {return _OutSubFolder;}
      
      std::string OutSubsubFolder() {return _OutSubsubFolder;} 
      void SetOutSubsubFolder(std::string string) { _OutSubsubFolder = string;}
      const std::string & OutSubsubFolder() const {return _OutSubsubFolder;}
      
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



