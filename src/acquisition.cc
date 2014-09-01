/**
 * \file
 * \brief Implementation of acquisition methods
 */

#include<string>
#include "acquisition.h"

namespace TCT {

  bool acquisition_single::Read(FILE *infile, uint32_t iFile){

    int debug = 8;
    // read input file
    if(debug>9) std::cout << "start reading"  << std::endl;

    int ret=0;
    float t_off=0;
    float Polarity = -1.;
    if(BiasVolt() < 0) Polarity = 1.;


    // read header to dummy
    Char_t dummy[20];
    for(Int_t k=0; k<30;k++){
      ret = fscanf(infile,"%s",&dummy);
      //if(debug>9) std::cout << k << " " << dummy << "    ";
    }

    // read pulse waveform
    float in_v, in_t;
    uint32_t counter = 0;
    while(1)
    {
      //std::cout << "\n4 " << std::endl;
      ret = fscanf(infile,"%g",&in_t);
      ret = fscanf(infile,"%g",&in_v);
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
      volt.push_back(Polarity*in_v);
      counter++;
    }

    SetNsamples(counter);


    // fill histogram and add to Object array // !! next 6 lines to be moved to the constructor
    Char_t buffername [50];
    Char_t buffername2 [50];
    sprintf (buffername, "Pulse_%d", (int)iFile);
    sprintf (buffername2, "Pulse_%d-10", (int)iFile);

    
    _H_acquisition = new TH1F(buffername, buffername, Nsamples(), time[0]-SampleInterval()*0.5, time[Nsamples()-1]+SampleInterval()*0.5);
    _H_acquisitionFILTERED = new TH1F(buffername2, buffername2, Nsamples(), time[0]-SampleInterval()*0.5, time[Nsamples()-1]+SampleInterval()*0.5);


    /* !! overload Read() method with another Read(...,TCT::acquisition avg) method, that writes this into avg object
    if(!hPulse) hPulse = new TH2F("acqs","acqs",Nsamples,         time[0]-SampleInterval()*0.5-50, time[Nsamples()-1]+SampleInterval()*0.5-50,560,-0.05,0.4);
    if(!hProf) hProf = new TProfile("avgpulses","avgpulses",Nsamples, pulse->t[0]-SampleInterval()*0.5-50, pulse->t[Nsamples-1]+SampleInterval()*0.5-50,-0.2,0.2," ");
    if(!hProffilt) hProffilt = new TProfile("avgpulses_filt","avgpulses-comb from filter",Nsamples, pulse->t[0]-SampleInterval()*0.5-50, pulse->t[Nsamples-1]+SampleInterval()*0.5-50,-0.2,0.4," ");

    if(!h_delay_width) h_delay_width = new TH2F("delay_width","delay_width",120,0,120,100,0,50);
    if(!h_ampl_width)  h_ampl_width  = new TH2F("ampl_width","ampl_width",100,0,0.4,100,0,50);
    if(!h_delay_ampl)  h_delay_ampl  = new TH2F("delay_ampl","delay_ampl",120,0,120,100,0,0.4);
    if(!h_rise1090_ampl)  h_rise1090_ampl  = new TH2F("rise1090_ampl","rise1090_ampl",100,0,20,100,0,0.4);
    if(!h_noise)	h_noise	      = new TH1F("noise","noise",100,0,0.01);
    */

    SetName("SingleAcq");
    //     pulse->Print();
    //AllTests->Add((BMTest *)pulse);	 	 


    //std::cout << "->return " << std::endl;
    return kTRUE;
  }

  void acquisition_single::SetName(std::string name){
    
    _name = name;

  }
}
