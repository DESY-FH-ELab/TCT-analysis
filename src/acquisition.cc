/**
 * \file
 * \brief Implementation of acquisition methods
 */

// STD includes
#include<string>

// TCT includes
#include "acquisition.h"
#include "analysis.h"
#include "util.h"

// ROOT includes
#include "TMath.h" 
#include "TRandom3.h" 

//#define DEBUG

namespace TCT {

  bool acquisition_single::Read(FILE *infile, uint32_t iFile){

#ifdef DEBUG 
    std::cout << "start ACQ_single::Read"  << std::endl;
#endif


    int ret=0;
    float t_off=0;

    Char_t dummy[20];
    // pre-parse to find polarity
    double in_v, in_t;
    float minvolt = 1.;
    float maxvolt = -1.;
    uint32_t counter = 0;
    for(Int_t k=0; k<30; k++){
      ret = fscanf(infile,"%s",dummy);
      //std::cout << "ret = " << ret << " k = " << k << " " << " dummy = " << dummy << "    ";
    }
    while(1)
    {
      ret = fscanf(infile,"%lf",&in_t);
      ret = fscanf(infile,"%lf",&in_v);
      if(feof(infile)) break;
      if(in_v > maxvolt) maxvolt = in_v;
      if(in_v < minvolt) minvolt = in_v;

      if (ret<=0) 
      {
	std::cout << "read error voltage block at position i = " << counter << "\n";
	exit(1);
      }
      counter++;

    }


    if(maxvolt > fabs(minvolt)) SetPolarity(1.);
    else SetPolarity(-1.);
    //if(BiasVolt() < 0) Polarity = 1.;

    rewind(infile);
    // read header to dummy again
    for(Int_t k=0; k<30; k++){
      ret = fscanf(infile,"%s",dummy);
      //if(debug>9) std::cout << k << " " << dummy << "    ";
    }

    // read pulse waveform
    counter = 0;
    while(1)
    {
      //std::cout << "\n4 " << std::endl;
      ret = fscanf(infile,"%lf",&in_t);
      ret = fscanf(infile,"%lf",&in_v);
      if(feof(infile)) break;

      //if(counter < 5) std::cout << counter << " in_t: " << in_t;
      //if(counter < 5) std::cout << "  " << counter << " in_v: " << in_v << " ";
      //if(counter < 5) std::cout << SampleInterval();

      if (ret<=0) 
      {
	std::cout << "read error voltage block at position i = " << counter << "\n";
	exit(1);
      }

      time.push_back(t_off+counter*SampleInterval());
      volt.push_back(Polarity()*in_v);
      counter++;
    }

    SetNsamples(counter);


    // fill histogram and add to Object array // !! next 6 lines to be moved to the constructor
    Char_t buffername [50];
    Char_t buffername2 [50];
    sprintf (buffername, "Pulse_%d", (int)iFile);
    sprintf (buffername2, "Pulse_%d-FILTERED", (int)iFile);


    _H_acquisition = new TH1F(buffername, buffername, Nsamples(), time[0]-SampleInterval()*0.5, time[Nsamples()-1]+SampleInterval()*0.5);
    _H_acquisitionFILTERED = new TH1F(buffername2, buffername2, Nsamples(), time[0]-SampleInterval()*0.5, time[Nsamples()-1]+SampleInterval()*0.5);


    this->SetName("SingleAcq");
    //this->PrintAcq();

#ifdef DEBUG 
    std::cout << "end ACQ_single::Read"  << std::endl;
#endif


    return kTRUE;
  }

  void acquisition_single::SetName(std::string name){

    _Name = name;

  }

  void acquisition_single::PrintAcq(){

    for( int i = 0; i < Nsamples(); i++) std::cout << " " << i << ": " << time[i] << " " << volt[i] << std::endl;

  }

  void acquisition_single::FillNtuple(TCT::acquisition_avg *acqAvg) {

#ifdef DEBUG 
    std::cout << "start ACQ_single::FillNtuple " << std::endl;
#endif

    float par[9];
    par[0] = _Rise;
    par[1] = _Rise1090;
    par[2] = _Fall;
    par[3] = _Width;
    par[4] = _Delay;
    par[5] = _Delayfilt;
    par[6] = _Maxamplitude;
    par[7] = _Avg;
    par[8] = _S2nval;

    (acqAvg->N_tuple())->Fill(par);

#ifdef DEBUG 
    std::cout << "end ACQ_single::FillNtuple " << std::endl;
#endif

    return;
  }

  void acquisition_single::GetOffsetNoise(uint32_t iAcq, TCT::acquisition_avg *acqAvg){

#ifdef DEBUG 
    std::cout << "start ACQ_single::GetOffsetNoise " << std::endl;
#endif

    //if(debug > 0) std::cout << "start GON" << std::endl;
    //std:: cout << Nsamples() << std::endl;
    //std:: cout << Nsamples_start() << std::endl;
    //std:: cout << Nsamples_end() << std::endl;

    float mean, mean_end, rms, rms_end;

    mean = .0;
    for (uint32_t i = 0; i < Nsamples_start(); i++) 
      mean += volt[i];
    mean /= (float)Nsamples_start();
    //std::cout << " mean = " << mean << std::endl;
    //if(Nsamples_start() > 0 && Nsamples_start() < Nsamples()) mean = TMath::Mean(Nsamples_start(), &volt[0]);

    rms = .0;
    for (uint32_t i = 0; i < Nsamples_start(); i++) 
      rms += (volt[i]-mean)*(volt[i]-mean);
    rms /= (float)Nsamples_start();
    rms = TMath::Power(rms,0.5);
    //if(Nsamples_start() > 0 && Nsamples_start() < Nsamples()) rms = TMath::RMS(Nsamples_start(), &volt[0]);

#ifdef DEBUG 
    std::cout << " Baseline offset = " << mean << std::endl;
    std::cout << " Baseline rms = " << rms << std::endl;
#endif

    mean_end = .0;
    for (uint32_t i = Nsamples() - Nsamples_end(); i < Nsamples(); i++) 
      mean_end += volt[i];
    mean_end /= (float)Nsamples_end();

    rms_end = .0;
    for (uint32_t i = Nsamples() - Nsamples_end(); i < Nsamples(); i++) 
      rms_end += (volt[i]-mean_end)*(volt[i]-mean_end);
    rms_end /= (float)Nsamples_end();
    rms_end = TMath::Power(rms_end,0.5);

#ifdef DEBUG 
    std::cout << " Baseline_end offset = " << mean_end << std::endl;
    std::cout << " Baseline_end rms = " << rms_end << std::endl;
#endif

    SetOffset(mean);
    SetNoise(rms);

    SetOffset_end(mean_end);
    SetNoise_end(rms_end);



    (acqAvg->H_noise())->Fill(Noise());

    (acqAvg->G_noise_evo())->SetPoint(iAcq,iAcq,Noise());
    // std::cout << "Noise: " << Noise << std::endl;

#ifdef DEBUG 
    std::cout << "end ACQ_single::GetOffsetNoise " << std::endl;
#endif
    return;

  }

  void acquisition_single::FillHacqs(){

#ifdef DEBUG 
    std::cout << "start ACQ_single::FillHacqs " << std::endl;
#endif

    for (Int_t j=0; j< Nsamples(); j++) (Hacq())->SetBinContent(j+1, volt[j] - Offset());
    for (Int_t j=0; j< Nsamples(); j++) {
      if( j <  HalfFilterwidth()) 
	(HacqFILTERED())->SetBinContent(j+1, (Hacq())->GetBinContent(j+1));
      if( j >= HalfFilterwidth() && j < Nsamples() - HalfFilterwidth()) 
	(HacqFILTERED())->SetBinContent(j+1, (Hacq())->Integral(j-HalfFilterwidth(),j+HalfFilterwidth())/(2.*HalfFilterwidth()));
      if( j >  Nsamples() - HalfFilterwidth()) 
	(HacqFILTERED())->SetBinContent(j+1, (Hacq())->GetBinContent(j+1));
    }

#ifdef DEBUG 
    std::cout << "start ACQ_single::FillHacqs " << std::endl;
#endif

    return;
  }

  void acquisition_single::Fill2DHistos(TCT::acquisition_avg *acqAvg){

#ifdef DEBUG 
    std::cout << "start ACQ_single::Fill2DHistos " << std::endl;
#endif

    (acqAvg->H2_delay_width())	->Fill(Delay(), Width());
    (acqAvg->H2_ampl_width())	->Fill(Maxamplitude(), Width());
    (acqAvg->H2_delay_ampl())	->Fill(Delay(), Maxamplitude());
    (acqAvg->H2_rise1090_ampl())	->Fill(Rise1090(), Maxamplitude());

#ifdef DEBUG 
    std::cout << "end ACQ_single::Fill2DHistos " << std::endl;
#endif

    return;
  }

  void acquisition_single::SignalFinder(TCT::acquisition_avg *acqAvg, float S2n_Cut, float Width_Cut, float Amplitude_Cut){

#ifdef DEBUG 
    std::cout << "start ACQ_single::SignalFinder " << std::endl;
    std::cout << " S2n_Cut = " << S2n_Cut << " Width_Cut = " << Width_Cut << " Amplitude_Cut = " << Amplitude_Cut << std::endl;
#endif

    Float_t s2n[Nsamples()];
    Int_t temp_start[1000];
    Int_t temp_end[1000];
    Float_t temp_sig[1000];
    Int_t start[1000];
    Int_t end[1000];
    Float_t sig[1000];
    Int_t AmpMaxPos[1000]; 
    Int_t temp_Found = 0;


    // init variables
    SetNFound(0);
    for (Int_t i=0; i<1000; i++) 
    {
      start[i] = -1;
      end[i] = -2;
      sig[i] = -999;
      AmpMaxPos[i] = -1;
    }
    for (Int_t i=0; i<1000; i++) 
    {
      temp_start[i] = -1;
      temp_end[i] = -2;
      temp_sig[i] = -999;
    }

    // calculate signal/noise ratio for each sample
    for (uint32_t i=0; i<Nsamples(); i++) 
      if (Noise() > 0) 
      {
	s2n[i] = volt[i] / Noise();
      } else {
	s2n[i] = 0;
	std::cout << " Noise <= 0 !!" << std::endl;
	return;
      }

    // find start and end of pulses
    Bool_t SigOn = kFALSE;
    Int_t Cnt =0;
    while (Cnt < Nsamples()-50)
    {
      if ((s2n[Cnt] + s2n[Cnt+1] + s2n[Cnt+2] + s2n[Cnt+3] + s2n[Cnt+4]) > 4.0*S2n_Cut &&  SigOn == kFALSE)
      {
	SigOn = kTRUE;
	temp_start[temp_Found] = Cnt;
	//std::cout << "TempStart: " << Cnt;
	temp_Found++;
      } else if ((s2n[Cnt] + s2n[Cnt+1] + s2n[Cnt+2] + s2n[Cnt+3] + s2n[Cnt+4] + s2n[Cnt+5]) < 5.0*S2n_Cut && SigOn == kTRUE) {
	temp_end[temp_Found-1] = Cnt;
	//std::cout << "TempStop: " << Cnt << std::endl;
	SigOn = kFALSE;
      } 
      if (temp_Found >= 1000) break;
      Cnt++;
    }
    //std::cout << "#temp found = " << temp_Found << endl;

    // check minimum pulse Length and calculate amplitude
    for (Int_t i=0; i<temp_Found; i++) {

      if ( (temp_end[i]-temp_start[i]) >= Width_Cut/SampleInterval()) {

	// get amplitude
	Float_t old_amp = -999;
	Int_t MaxSigPos = -1;
	for (Int_t j=temp_start[i]; j<temp_end[i]; j++) {
	  if (volt[j]>old_amp) {
	    MaxSigPos = j; 
	    old_amp = volt[j];
	  }
	}
	temp_sig[i] = (volt[MaxSigPos-1]+volt[MaxSigPos]+volt[MaxSigPos+1])/3.;

	// check is signal is larger than cut
	if (temp_sig[i] > Amplitude_Cut) {
	  start[NFound()] = temp_start[i]-1;
	  end[NFound()] = temp_end[i];
	  sig[NFound()] = temp_sig[i];
	  AmpMaxPos[NFound()] = MaxSigPos;
	  SetNFound(NFound()+1);
	  if (NFound() > 1000) {
	    std::cout << "   *** FOUND 1000 temporary pulses, breaking! " << std::endl;
	    break;
	  }
	}

      }
    }

    // check for biggest signal among all
    Float_t old_amp = -999;
    Int_t MaxSigLoc = -1;
    Float_t charge = 0.0;
    for(Int_t i =0; i<NFound(); i++) {
      charge = (time[end[i]]-time[start[i]])*sig[i]*0.75;
      if (charge>old_amp) {
	MaxSigLoc = i;
	old_amp = charge;
      }   
    }


    if (NFound()) {

      SetDelay(time[start[MaxSigLoc]]);


      int count = start[MaxSigLoc]-20;
      while(HacqFILTERED()->GetBinContent(count) + HacqFILTERED()->GetBinContent(count+1) + HacqFILTERED()->GetBinContent(count+2) < 3.*S2n_Cut*Noise()){
	count++;
	if(count >  end[MaxSigLoc]) {
	  std::cout <<"problem" << std::endl;
	  break;
	}

      }
      SetDelayfilt(((float)count)/10.);

      SetWidth(time[end[MaxSigLoc]]-time[start[MaxSigLoc]]);


      SetMaxamplitude(sig[MaxSigLoc]);

      float avg = 0.0;
      count = 0;
      for(int i = start[MaxSigLoc]+10; i < end[MaxSigLoc]-10; i++){
	if(volt[i]> 0.3*Maxamplitude()){
	  avg += volt[i];
	  count++;
	}

      }
      SetAvg(avg/((float)count));

      float avgshort = 0;
      count = 0;
      for(int i = AmpMaxPos[MaxSigLoc]-5; i < AmpMaxPos[MaxSigLoc]+5; i++){
	//std::cout << i << " " << volt[i] << "   ";
	avgshort += volt[i];
	count++;
      }
      //std::cout << " sum: " << avgshort << std::endl;
      //std::cout << " count: " << count << std::endl;

      SetAvgshort(avgshort/((float)count));

      SetS2nval(Avg()/Noise());


      (acqAvg->G_s2n_evo())->SetPoint(iAcq(),iAcq(),S2nval());

      SetRise(time[AmpMaxPos[MaxSigLoc]]-time[start[MaxSigLoc]]);

      SetFall(time[end[MaxSigLoc]] - time[AmpMaxPos[MaxSigLoc]]);

      //calculate risetime

      int rise10 = -1;
      int rise90 = -1;

      bool found10 = false;
      bool found90 = false;
      float frac = 0.25;
      if(IsMIP()) frac = 0.49;
      for(int i = start[MaxSigLoc]; i < AmpMaxPos[MaxSigLoc]; i++){ // MIPs
	//std::cout << "prof->GetBinContent(i): " <<  i << "  " << prof->GetBinContent(i) << std::endl;

	if((volt[i] + volt[i+1] + volt[i+2] + volt[i+3]) > 4.*0.1*Avgshort() && !found10){
	  rise10 = i;
	  found10 = true;
	  //std::cout << " rise10: " << rise10;
	}
	if(found10 && (volt[i] + volt[i+1] + volt[i+2]) > 3.*0.9*Avgshort() && !found90){
	  rise90 = i;
	  found90 = true;
	  //std::cout << "  rise90: " << rise90;
	}
      }
      if(!found90) rise90 = AmpMaxPos[MaxSigLoc];
      if (((float)(rise90 - rise10))/10. > 0.) SetRise1090(((float)(rise90 - rise10))/10.);

      if(Rise1090() < 0.0){ 
	std::cout << " NEG RISE TIME !!! " << std::endl;
	std::cout << "rise10 = " <<  rise10 << std::endl;
	std::cout << "rise90 = " <<  rise90 << std::endl;
	std::cout << "rise1090 = " <<  Rise1090() << std::endl;
	std::cout << "Avgshort = " << Avgshort() << std::endl;
      }


    } else {
      std::cout << " No pulse found " << std::endl;
    }

    for (Int_t i=0; i<Nsamples()-5; i++){
      if (i>2 && i < (int)(Delay()/SampleInterval()) ) if (volt[i] < AmplNegEarly() ) SetAmplNegEarly((volt[i]+volt[i+1]+volt[i-1])/3.);
      if (i>2 && i < (int)(Delay()/SampleInterval()) ) if (volt[i] > AmplPosEarly() ) SetAmplPosEarly((volt[i]+volt[i+1]+volt[i-1])/3.);
      if ( i >  (int)((Delay() + 2.*Width() + 2.)/SampleInterval())) if (volt[i] < AmplNegLate() ) SetAmplNegLate((volt[i]+volt[i+1]+volt[i-1]+volt[i+2]+volt[i-2])/5.);
      if ( i >  (int)((Delay() + 2.*Width() + 2.)/SampleInterval())) if (volt[i] > AmplPosLate() ) SetAmplPosLate((volt[i]+volt[i+1]+volt[i-1]+volt[i+2]+volt[i-2]-5.*Offset())/5.);
      //std::cout << s2n << "	";
    } 

    TCT::util util;

    if(this->NFound() > 1) {
      std::cout << "more than one pulse found! " << std::endl;
      std::cout << *this << std::endl;
    }

#ifdef DEBUG 
    std::cout << "end ACQ_single::SignalFinder " << std::endl;
#endif


    return;
  } // end signal finder

  void acquisition_single::NoiseAdder(float additional_noise_rms){

    TRandom *r3 = new TRandom3();
    r3->SetSeed(0);
    double x;
    for(int i = 0; i < Nsamples(); i++){
      x = r3->Gaus(0,additional_noise_rms);
      volt[i] += x;
    }

    //repeat Noise calculation for smeared stuff
    float mean, rms;

    mean = .0;
    for (uint32_t i = 0; i < Nsamples_start(); i++) 
      mean += volt[i];
    mean /= (float)Nsamples_start();

    rms = .0;
    for (uint32_t i = 0; i < Nsamples_start(); i++) 
      rms += (volt[i]-mean)*(volt[i]-mean);
    rms /= (float)Nsamples_start();
    rms = TMath::Power(rms,0.5);

    SetOffset(mean);
    SetNoise(rms);

    return;
  }

  void acquisition_single::JitterAdder(float additional_jitter_rms){

    TRandom *r3 = new TRandom3();
    r3->SetSeed(0);
    double t;
    t = r3->Gaus(0,additional_jitter_rms);

    SetDelay(Delay() + t);
    return;
  }
}
