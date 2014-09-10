/**
 * \file
 * \brief Implementation of analysis methods
 */

// STD includes
#include<string>

// TCT includes
#include "analysis.h"
#include "acquisition.h"
//#include "util.h"

//  ROOT includes
#include "TSystem.h"
#include "TFile.h"

namespace TCT {

  void analysis::SetParameters(std::map<std::string, std::string> id_val){

    #ifdef DEBUG 
    std::cout << " start ANA::SetParameterd: Read and set parameters from input map" << std::endl; 
    #endif

    for( auto i : id_val){
      //if(i.first == "") _ = atof((i.second).c_str());
      if(i.first == "AmplNegLate_Cut") _AmplNegLate_Cut = atof((i.second).c_str());
      if(i.first == "AmplPosLate_Cut") _AmplPosLate_Cut = atof((i.second).c_str());
      if(i.first == "AmplNegEarly_Cut") _AmplNegEarly_Cut = atof((i.second).c_str());
      if(i.first == "AmplPosEarly_Cut") _AmplPosEarly_Cut = atof((i.second).c_str());
      if(i.first == "Outfolder") _OutFolder = i.second;
      if(i.first == "Noise_Cut") _Noise_Cut = atof((i.second).c_str());
      if(i.first == "NoiseEnd_Cut") _NoiseEnd_Cut = atof((i.second).c_str());
      if(i.first == "S2n_Cut") _S2n_Cut = atof((i.second).c_str());
      if(i.first == "S2n_Ref") _S2n_Ref = atof((i.second).c_str());
      if(i.first == "PrintEvent") _PrintEvent = atoi((i.second).c_str());
    }

    #ifdef DEBUG 
    std::cout << " end ANA::SetParameters" << std::endl; 
    #endif

    return;
  }

  bool analysis::AcqsSelecter(TCT::acquisition_single *acq){

    #ifdef DEBUG
    std::cout << " start ANA::AcqsSelecter" << std::endl;
    #endif

    bool ok = true;

    if(acq->SelectionRan()) return acq->Select();
    else {
      if (acq->Maxamplitude() < .0) {
	ok = kFALSE;
	std::cout << "no pulse found " << std::endl;
	return ok;
      }
      if(acq->Noise() > Noise_Cut()) {
	ok = false;
	std::cout << " acq too noisy: " << acq->Noise() << " Noise cut = " << Noise_Cut() << std::endl;
      }
      if (acq->Noise_end() > NoiseEnd_Cut()) {
	ok = kFALSE;
	std::cout << "pulse end too noisy -> pick-up or acquisition window too narrow" << std::endl;
	return ok;
      }
      if (acq->S2nval() < S2n_Cut()){
	ok = kFALSE;
	std::cout << "s2n too small: " << acq->S2nval() << std::endl;
	return ok;
      }
      //if (acq->Width() < Width_Cut()){
      //  ok = kFALSE;
      //  std::cout << "width too small: " << acq->Width()" << std::endl;
      //  return ok;
      //} // check width already during SignalFinder
      if (acq->AmplNegLate() < AmplNegLate_Cut()) {
	ok = kFALSE;
	std::cout << "pulse " << acq->iAcq() << " had big neg component after pulse" << acq->AmplNegLate() <<  std::endl;
	return ok;
      }
      if (acq->AmplPosLate() > AmplPosLate_Cut()) {
	ok = kFALSE;
	std::cout << "pulse " << acq->iAcq() << " had big pos component after pulse" << acq->AmplPosLate() <<  std::endl;
	return ok;
      }
      if (acq->AmplPosEarly() > AmplPosEarly_Cut()) {
	ok = kFALSE;
	std::cout << "pulse " << acq->iAcq() << " had big pos component before pulse" << acq->AmplPosEarly() <<  std::endl;
	return ok;
      }
      if (acq->AmplNegEarly() < AmplNegEarly_Cut()) {
	ok = kFALSE;
	std::cout << "pulse " << acq->iAcq() << " had big neg component before pulse" << acq->AmplNegEarly() <<   std::endl;
	return ok;
      }
    }

    acq->SetSelectionRan(true);
    acq->SetSelect(ok);

    return ok;
  }

  void analysis::AcqsAnalyser(TCT::acquisition_single *acq, uint32_t iAcq, TCT::acquisition_avg *acqAvg){

    #ifdef DEBUG
    std::cout << "start ANA::AcqsAnalyser" << std::endl;
    #endif

    //std::cout << "GOF" << std::endl;
    acq->GetOffsetNoise(iAcq, acqAvg);

    //std::cout << "FHacqs" << std::endl;
    acq->FillHacqs();

    //std::cout << "SF" << std::endl;
    acq->SignalFinder(acqAvg, S2n_Cut(), Width_Cut(), Amplitude_Cut());

    if(acq->iAcq() == PrintEvent()) std::cout << *acq << std::endl;

    //std::cout << "FillNt" << std::endl;
    acq->FillNtuple(acqAvg); // !! this is w/o selection cuts

    //std::cout << "Fill2DHs" << std::endl;
    acq->Fill2DHistos(acqAvg); // !! this is w/o selection cuts


    // amount of delayfilt depends on steepness, hence needs shifting now

    /*float meanfilt = 0.;
      for(Int_t i = 0; i<300; i++){
    //cout << ",  " << hProffilt->GetBinContent(i);
    meanfilt += hProffilt->GetBinContent(i+1);
    }
    meanfilt = meanfilt / 300.;
    cout << "Meanfilt = " << meanfilt << endl;

    float rmsfilt = 0.;
    for(Int_t i = 0; i<300; i++) rmsfilt += (hProffilt->GetBinContent(i+1) - meanfilt)*(hProffilt->GetBinContent(i+1) - meanfilt);
    rmsfilt /= 300.;
    cout << "rmsfilt = " << rmsfilt << endl;
    float noisefilt = TMath::Power(rmsfilt,0.5);
    cout << "noisefilt = " << noisefilt << endl;

    int count = 300;
    float limit = 3.*2*sqrt(30)*noisefilt + meanfilt;
    cout << "limit = " << limit << endl;
    while(hProffilt->GetBinContent(count) + hProffilt->GetBinContent(count+1) + hProffilt->GetBinContent(count+2) < limit){
    count++;
    if(count >  1000) {
    cout <<"problem for hProffilt" <<endl;
    break;
    }

    }
    cout << "count = " << count << endl;
    float delayfilt2 = ((float)count-500.)/10.;
    cout << "delayfilt2 = " << delayfilt2 << endl;



    */   

    #ifdef DEBUG
    std::cout << "end ANA::AcqsAnalyser" << std::endl;
    #endif

    return;	
  }
 
  void analysis::AcqsProfileFiller(TCT::acquisition_single *acq, TCT::acquisition_avg *acqAvg) {

    #ifdef DEBUG
    std::cout << "start ANA::AcqsProfileFiller" << std::endl;
    #endif

    if (acq->Select()){
      Float_t tmp_t = -1.0;
      Float_t tmp_tfilt = -1.0;
      Float_t tmp_v = -1.0;

      for(Int_t j = 0; j < acq->Nsamples(); j++) {
	tmp_t = acq->time[j] - acq->Delay();
	tmp_tfilt = acq->time[j] - acq->Delayfilt();
	tmp_v = acq->volt[j] - acq->Offset();

	acqAvg->H2_acqs2D()->Fill(tmp_t,tmp_v,1.);
	acqAvg->Profile()->Fill(tmp_t,tmp_v,1.);
	acqAvg->ProfileFILTERED()->Fill(tmp_tfilt,tmp_v,1.); // !! ProfileFILTERED is not filtered, but shifted by certain amount with calculation based on filtered signal, for low SNR only
      }

    } 

    #ifdef DEBUG
    std::cout << "end ANA::AcqsProfileFiller" << std::endl;
    #endif

    return;

  }

    void analysis::AcqsWriter(TCT::sample *sample, std::vector<TCT::acquisition_single> *allAcqs, TCT::acquisition_avg *acqAvg){

    #ifdef DEBUG
    std::cout << "start ANA::AcqsWriter" << std::endl;
    #endif

    std::string outfolder	= OutFolder(); 
    std::string outsample	= sample->SampleID();
    std::string outtemp		= std::to_string((int)acqAvg->Temp());
    std::string outpolarity = "+";
    if(acqAvg->Polarity() > .0) outpolarity = "-";
    std::string outvolt 	= std::to_string((int)acqAvg->BiasVolt());


    std::string outpath  = outfolder + "/" + outsample + "/" + outtemp + "K";
    std::string outpath1 = outfolder + "/" + outsample + "/";
    std::string pathandfilename = outpath  + "/" + outsample + "_" + outtemp + "K_" + outpolarity + outvolt + "V.root";

    gSystem->MakeDirectory(outpath1.c_str());
    //if(gSystem->MakeDirectory(outpath.c_str()) == -1) std::cout << "couldnt create directory" << std::endl;
    gSystem->MakeDirectory(outpath.c_str());

    std::cout << "\n *** outfile written to: " << pathandfilename << " *** " << std::endl;

    TFile* f_rootfile = new TFile(pathandfilename.c_str(),"RECREATE","TCTanalyser");

    f_rootfile->cd();
    acqAvg->N_tuple()->Write();
    acqAvg->H2_acqs2D()->Write();
    acqAvg->Profile()->Write();
    //acqAvg->ProfileFILTERED()->Write();
    acqAvg->H2_delay_width()->Write();
    acqAvg->H2_ampl_width()->Write();
    acqAvg->H2_delay_ampl()->Write();
    acqAvg->H2_rise1090_ampl()->Write();
    acqAvg->H_noise()->Write();
    acqAvg->G_noise_evo()->Write();
    acqAvg->G_s2n_evo()->Write();

    // now write the single pulses into the file !! add option -save_all for switch
    f_rootfile->mkdir("single_acqs");
    f_rootfile->cd("single_acqs");

    //std::cout << allAcqs->at(0) << std::endl;
    for(uint32_t i_acq = 0; i_acq < allAcqs->size(); i_acq++){

      TCT::acquisition_single* acq = &allAcqs->at(i_acq);
      acq->Hacq()->Write();

      //std::cout << *acq << std::endl;


    }
    //for ( TCT::acquisition_single acq : &allAcqs) {
    //  std::cout << acq << std::endl;
    //}

    f_rootfile->Close();

    #ifdef DEBUG
    std::cout << "end ANA::AcqsWriter" << std::endl;
    #endif

    return;

  }

}
