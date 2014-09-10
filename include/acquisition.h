/**
 * \file
 * \brief Definition of the sample class.
 */

#ifndef __ACQUISITION_H__
#define __ACQUISITION_H__ 1

// STD includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TGraph.h"

// TCT includes
//#include "analysis.h"

namespace TCT {

  // \brief Abstract class for acquisition, acq_single and  acq_multi inherit from it

  class acquisition_base {

    private :

      bool _IsEdgeTCT;
      double _x, _y; // do we need local and global coordinates ?? in um
      double _SNR;
      double _AmpRefLaser;
      bool _IsManyPulseStructure;
      bool _IsMIP;
      uint32_t _Nsamples;
      uint32_t _Nsamples_start;
      uint32_t _Nsamples_end;

    protected:

      float _Temp;	// in K
      float _BiasVolt;	// in Volt
      float _SampleInterval;	// in ns
      float _PrePulseInterval;	// in ns. !! By how much should be translated can be inferred from Delay(), or from the trigger position written in the scope header
      float _Polarity;

    public :
    
      acquisition_base() :
        _Temp(295.),
	_BiasVolt(500.),
	_Nsamples(-1),
	_Nsamples_start(50),
	_Nsamples_end(50),
	_SampleInterval(0.1),
	_PrePulseInterval(20.),
	_Polarity(.0)
      {};

      /*acquisition_base(float bias) :
	_BiasVolt(bias),
	_Nsamples(-1),
	_Nsamples_start(10),
	_Nsamples_end(10),
	_SampleInterval(0.1)
	{};*/

      acquisition_base(uint32_t nsamples) :
        acquisition_base() 
      {   
	_IsMIP = false; // !! Get this from measurement class
	_Nsamples = nsamples;
      };

      // Default copy constructer should be fine
      acquisition_base(const acquisition_base &)               = default;
      acquisition_base & operator = (const acquisition_base &) = default;

      // Dectructor
      ~acquisition_base() = default;


      virtual double GetSNR() const = 0;

      void SetNsamples(uint32_t nsamples) { _Nsamples = nsamples;}
      uint32_t Nsamples(){ return _Nsamples;}
      const uint32_t & Nsamples() const{ return _Nsamples;}

      void SetNsamples_start(uint32_t nsamples) { _Nsamples_start = nsamples;}
      uint32_t Nsamples_start(){ return _Nsamples_start;}
      const uint32_t & Nsamples_start() const{ return _Nsamples_start;}

      void SetNsamples_end(uint32_t nsamples) { _Nsamples_end = nsamples;}
      uint32_t Nsamples_end(){ return _Nsamples_end;}
      const uint32_t & Nsamples_end() const{ return _Nsamples_end;}

      float SampleInterval(){ return _SampleInterval;}
      void SetSampleInterval(float interval){ _SampleInterval = interval;}
      const float & SampleInterval() const{ return _SampleInterval;}

      float PrePulseInterval(){ return _PrePulseInterval;}
      void SetPrePulseInterval(float interval){ _PrePulseInterval = interval;}
      const float & PrePulseInterval() const{ return _PrePulseInterval;}

      std::vector<double> volt; // ?? encapsulate?
      std::vector<double> time;

      void SetBiasVolt(float volt) { _BiasVolt = volt;}
      float BiasVolt() {return _BiasVolt;}
      const float & BiasVolt() const { return _BiasVolt;}

      void SetTemp(float temp) { _Temp = temp;} 
      float Temp() {return _Temp;} 
      const float & Temp() const { return _Temp;}

      bool IsMIP() { return _IsMIP;}
      void SetIsMIP(bool Is) { _IsMIP = Is;}

      void SetPolarity(float pol) { _Polarity = pol;}
      float Polarity() {return _Polarity;}
      const float & Polarity() const { return _Polarity;}

  }; // end of class acquisition_base


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class acquisition_avg : public acquisition_base {

    private : 
      TProfile*	_Profile;
      TProfile*	_ProfileFILTERED;
      TH2F*	_H2_acqs2D;
      TH2F*	_H2_delay_width;
      TH2F*	_H2_ampl_width;
      TH2F*	_H2_delay_ampl;
      TH2F*	_H2_rise1090_ampl;
      TNtuple*	_N_tuple;
      TH1F*	_H_noise;
      TGraph*	_G_noise_evo;
      TGraph*	_G_s2n_evo;


    public :

      //acquisition_avg() = default;

      acquisition_avg(uint32_t nsamples) :
	acquisition_base(nsamples){
	  _H_noise 	= new TH1F("noise","noise",100,0,0.01);
	  _G_noise_evo 	= new TGraph(1000);
	  _G_noise_evo->SetNameTitle("Noise Evolution","Noise Evolution");
	  _G_noise_evo->SetMarkerStyle(2);
	  _G_s2n_evo 	= new TGraph(1000);
	  _G_s2n_evo->SetNameTitle("S2Noise Evolution","S2Noise Evolution");
	  _G_s2n_evo->SetMarkerStyle(2);
	  _N_tuple 	= new TNtuple("ntuple","pulse ntuple","rise:rise1090:fall:width:delay:delayfilt:ampl:avg:s2nval");
	  _H2_acqs2D 	= new TH2F("acqs2D","acqs2D", Nsamples(), -SampleInterval()*0.5 - PrePulseInterval(), (Nsamples()-1-0.5)*SampleInterval() - PrePulseInterval(), 5000, -1., 1.); 
	  _Profile 	= new TProfile("avgAcq","avgAcq",Nsamples(), -SampleInterval()*0.5 - PrePulseInterval(), (Nsamples()-1-0.5)*SampleInterval() - PrePulseInterval() ,-1,1," ");
	  _ProfileFILTERED = new TProfile("avgAcq_f","avgAcq_f",Nsamples(), -SampleInterval()*0.5-PrePulseInterval(), (Nsamples()-1-0.5)*SampleInterval(),-1,1," ");
	  _H2_delay_width = new TH2F("delay_width","delay_width",120,0,120,100,0,50);
	  _H2_ampl_width  = new TH2F("ampl_width","ampl_width",100,0,0.4,100,0,50);
	  _H2_delay_ampl  = new TH2F("delay_ampl","delay_ampl",120,0,120,100,0,0.4);
	  _H2_rise1090_ampl  = new TH2F("rise1090_ampl","rise1090_ampl",100,0,20,100,0,0.4);
	}

      //_Nsamples(1996),
      //{}

      // Default copy constructer should be fine
      acquisition_avg(const acquisition_avg &)               = default;
      acquisition_avg & operator = (const acquisition_avg &) = default;

      // Dectructor
      ~acquisition_avg() = default;

      virtual double GetSNR() const {
	return 1.0;
      }

      TProfile*	Profile() 	{ return _Profile;}
      TProfile*	ProfileFILTERED()	{ return _ProfileFILTERED;}
      TH2F*	H2_acqs2D()	{ return _H2_acqs2D;}
      TH2F*	H2_delay_width()	{ return _H2_delay_width;}
      TH2F*	H2_ampl_width()	{ return _H2_ampl_width;}
      TH2F*	H2_delay_ampl()	{ return _H2_delay_ampl;}
      TH2F*	H2_rise1090_ampl()	{ return _H2_rise1090_ampl;}
      //TH2F*	H2__() { return _H2__;}
      TNtuple*	N_tuple()	{ return _N_tuple;}
      TH1F*	H_noise()	{ return _H_noise;}
      TGraph*	G_noise_evo()	{ return _G_noise_evo;}
      TGraph*	G_s2n_evo()	{ return _G_s2n_evo;}

  }; // end of acquisition_avg implementation


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // \brief Acquisition with single pulse implementation

  class acquisition_single : public acquisition_base {

    private :

      std::string _Name;
      uint32_t _iAcq;
      float _Maxamplitude;	// this is the amplitude of a single pulse, at MaxSigLocation
      float _AmplNegEarly;		// most neg value before rising edge
      float _AmplPosEarly;		// most pos value before rising edge
      float _AmplNegLate;	// most neg value after falling edge
      float _AmplPosLate;	// most neg value after falling edge
      float _Avg;		// between edges for alpha 
      float _Avgshort;		// avg shortly after risind edge (beta,p,...)
      float _Delay;		// from acq start to start of RE, start of RE is defined by S2nRef in param.h
      float _Delayfilt;		// from acq start to start of RE for filtered pulses
      float _Width;		// from 50% of RE to 50% of FE
      float _Rise;		// from start to position of maxamplitude 
      float _Rise1090;		// rise time from 10% to 90%
      float _Fall;		// from position of maxamp to end of pulse, end of FE is defined by S2nRef in param.h
      float _S2nval;		// avg / Noise
      float _Offset;		// baseline offset before RE
      float _Offset_end;	// baseline offset after FE
      float _Noise;		// Noise before RE
      float _Noise_end; 	// Noise after FE
      TH1F *_H_acquisition;		// TH1F to store the data
      TH1F *_H_acquisitionFILTERED;	// TH1F to store the data, filtered

      int _HalfFilterwidth;

      bool _SelectionRan;
      bool _Selected;

      uint32_t _NFound;


    public :
      acquisition_single() : 
	_Name("single"),
	acquisition_base(500.),
	_SelectionRan(false),
	_NFound(0){};

      acquisition_single(uint32_t iAcq) : 
	acquisition_base(500.),
	_Name("single"),
	_iAcq(iAcq),
	_Maxamplitude(-1.),
	_AmplNegEarly(1.),
	_AmplPosEarly(-1.),
	_AmplNegLate(1.),
	_AmplPosLate(-1.),
	_Avg(-1.),
	_Avgshort(-1.),
	_Delay(-1.),
	_Delayfilt(-1.),
	_Width(-1.),
	_Rise(-1.),
	_Rise1090(-1.),
	_Fall(-1.),
	_SelectionRan(false),
	_NFound(0){};


      // Default copy constructer should be fine
      acquisition_single(const acquisition_single &)               = default;
      acquisition_single & operator = (const acquisition_single &) = default;

      // Dectructor
      ~acquisition_single() = default;


      virtual double GetSNR() const {
	return 1.0;
      }

      void Print();
      void GetOffsetNoise(uint32_t iEv, TCT::acquisition_avg *avg);
      bool Select() { return _Selected;}
      bool SelectionRan() { return _SelectionRan;}
      void SetSelectionRan(bool select) { _SelectionRan = select;}
      void SetSelect(bool select) { _Selected = select;}

      void DrawPulse();
      bool Read(FILE *infile, uint32_t iFile);
      bool Read(FILE *infile, uint32_t iFile, TCT::acquisition_avg *avg); 
      void FillNtuple(TCT::acquisition_avg *avg);
      void SignalFinder(TCT::acquisition_avg *avg, float, float, float );
      void SignalManipulator();
      void ClearStruct();
      void SetName(std::string name);
      std::string & Name() { return _Name;}
      const std::string & Name() const { return _Name;}

      void SetiAcq(uint32_t iacq) { _iAcq = iacq;}
      uint32_t iAcq() { return _iAcq;}
      const uint32_t & iAcq() const { return _iAcq;}

      TH1F* Hacq() { return _H_acquisition;}
      TH1F* HacqFILTERED() { return _H_acquisitionFILTERED;}

      float Offset() { return _Offset;}
      void SetOffset(float offset) { _Offset = offset;}
      const float & Offset() const { return _Offset;}

      float Offset_end() { return _Offset_end;}
      void SetOffset_end(float offset) { _Offset_end = offset;}
      const float & Offset_end() const { return _Offset_end;}

      float Noise() { return _Noise;}
      void SetNoise(float noise) { _Noise = noise;}
      const float & Noise() const { return _Noise;}

      float Noise_end() { return _Noise_end;}
      void SetNoise_end(float noise) { _Noise_end = noise;}
      const float & Noise_end() const { return _Noise_end;}

      void SetHalfFilterwidth(uint32_t halfwidth){ _HalfFilterwidth = halfwidth;}
      uint32_t HalfFilterwidth(){ return _HalfFilterwidth;}

      void FillHacqs();
      void Fill2DHistos(TCT::acquisition_avg *avg);

      float Delay() { return _Delay;}
      void SetDelay(float delay) { _Delay = delay;}
      const float & Delay() const { return _Delay;}

      float Delayfilt() { return _Delayfilt;}
      void SetDelayfilt(float delay) { _Delayfilt = delay;}

      float Width() { return _Width;}
      void SetWidth(float width) { _Width = width;}
      const float & Width() const { return _Width;}

      float Maxamplitude() { return _Maxamplitude;}
      void SetMaxamplitude(float amp) { _Maxamplitude = amp;}
      const float & Maxamplitude() const { return _Maxamplitude;}

      float AmplNegEarly() { return _AmplNegEarly;}
      void SetAmplNegEarly(float amp) { _AmplNegEarly = amp;}

      float AmplPosEarly() { return _AmplPosEarly;}
      void SetAmplPosEarly(float amp) { _AmplPosEarly = amp;}

      float AmplNegLate() { return _AmplNegLate;}
      void SetAmplNegLate(float amp) { _AmplNegLate = amp;}

      float AmplPosLate() { return _AmplPosLate;}
      void SetAmplPosLate(float amp) { _AmplPosLate = amp;}

      float Rise() { return _Rise;}
      void SetRise(float time) { _Rise = time;}

      float Rise1090() { return _Rise1090;}
      void SetRise1090(float time) { _Rise1090 = time;}
      const float & Rise1090() const { return _Rise1090;}

      float Fall() { return _Fall;}
      void SetFall(float time) { _Fall = time;}
      const float & Fall() const { return _Fall;}

      float S2nval() { return _S2nval;}
      void SetS2nval(float time) { _S2nval = time;}
      const float & S2nval() const { return _S2nval;}

      float Avg() { return _Avg;}
      void SetAvg(float avg) { _Avg = avg;}

      float Avgshort() { return _Avg;}
      void SetAvgshort(float avg) { _Avgshort = avg;}

      uint32_t NFound() { return _NFound;}
      void SetNFound(uint32_t n) { _NFound = n;}
      const uint32_t & NFound() const { return _NFound;}

  }; // end of acquisition_single implementation

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class acquisition_multi_single : public acquisition_base {

    public :
      acquisition_multi_single() = default;//:
      //_Nsamples(1996),
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
      //_Nsamples(1996),
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

inline std::ostream & operator << (std::ostream & os, const TCT::acquisition_single & acq) {
  return os	<< "   This acq is " << acq.Name() << " #" << acq.iAcq()
    << "\n MaxAmpl " << acq.Maxamplitude()
    << "   Delay " << acq.Delay()
    << "   Width " << acq.Width()
    << "   Rise1090 " << acq.Rise1090()
    << "   Fall " << acq.Fall()
    << "   S2N ratio " << acq.S2nval()
    << "   Noise " << acq.Noise()
    << "   Noise_End " << acq.Noise_end()
    << "   Offset " << acq.Offset()
    << "   Offset_End " << acq.Offset_end()
    << "   NFound " << acq.NFound()
    //<< "    " << acq.()
    << std::endl;
}
#endif
