/**
 * \file
 * \brief Main program code.
 */

//  includes from standard libraries
#include <iostream>
#include <regex>
#include <fstream>

//  includes from TCT classes
#include "sample.h"
#include "config.h"
//#include "param.h"
#include "util.h"
#include "acquisition.h"
#include "measurement.h"
#include "analysis.h"

//  includes from ROOT libraries
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

//#define DEBUG

int main(int argc, char* argv[])
{
  std::cout << "\n  This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "\n" << std::endl;

  /*
  //TApplication theApp("App", 0, 0);
  TCanvas* c;  
  c = new TCanvas("Diamond Beam Monitor","Diamond Beam Monitor",1200,800);
  c->cd();
  gROOT->SetStyle("Plain");
  TGraph* g = new TGraph(3);
  g->SetPoint(0,1.,1.);
  g->SetPoint(1,2.,2.);
  g->SetPoint(2,3.,3.);
  g->Draw("apl");

  c->Modified();
  c->Update();
   */

  if(argc == 1){
    std::cout << " No root folder specified. Please execute with path to projcet: > ./tct-analysis -r /home/<user>/<my-path>/TCT-analysis/" << std::endl;
    return 1;
  }

  TCT::util util;

  std::string proj_folder = "default";
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-r")) {
      proj_folder = argv[++i];
      std::cout << "Root folder of TCT project is " << proj_folder << std::endl;
    }
    if (!strcmp(argv[i],"-af")) {
      std::ifstream ana_file (argv[++i]);
      std::cout << "Analysis file is " << argv[i] << std::endl;
      util.parse(ana_file);
    }
    // !! add check for certain vital options, if not passed, break

    // Maximum events:
    /*if (!strcmp(argv[i],"-e")) {
      max_events = atoi(argv[++i]);
      std::cout << "Decoding a maximum number of " << max_events << "
      events." << std::endl;
      }*/
  }

  uint32_t MaxAcqs= -1;
  for( auto i : util.ID_val()){
    if(i.first == "ProjectFolder") {
      if(proj_folder != "default") std::cout << " Project folder from command line overwritten! " << std:: endl;
      proj_folder = i.second;
    }
    if(i.first == "MaxAcqs") {
      MaxAcqs = atoi(i.second.c_str());
    }
  }

  if(proj_folder == "default") {
    std::cout << " project folder not specified, neither in command line nor in analysis file!\n\n   ***STOPPING" << std::endl;
    return 1;
    
  }

  TCT::analysis ana(util.ID_val());
  //std::cout << ana << std::endl;



  std::string DataFolder	= proj_folder + "/testdata/S57/295K/500V/";
  //std::string OutFolder	= proj_folder + "/results";
  std::string SensorFolder	= proj_folder + "/testSensor";

  TCT::sample dummyDUT2(SensorFolder);       // define DUT
  //std::cout << dummyDUT2 << std::endl;	// print basic parameters of the sample
  dummyDUT2.ReadSampleCard();	// read SampleCard and set parameters accordingly

  std::string sampleID = "S57";
  dummyDUT2.SetSampleID(sampleID); 
  //std::cout << dummyDUT2 << std::endl;

  std::vector<TCT::acquisition_single> AllAcqs;


  //TCT::param param;

  TCT::measurement meas(DataFolder);
  //std::cout << meas << std::endl;

  if(!meas.AcqsLoader(&AllAcqs, MaxAcqs)) return 1;//, 4); // change to take parameter from param
  // !! analyser functions belong to measurement class a.t.m., disentangle ?

  // now create instance of avg acquisition using Nsamples from loaded files
  TCT::acquisition_avg AcqAvg(AllAcqs[0].Nsamples());

  //now analyse all acquisitions
  int Nselected = 0;

  #ifdef DEBUG 
  std::cout << "Size of AllAcqs = " << AllAcqs.size() << std::endl;
  #endif

  for(uint32_t i_acq = 0; i_acq < AllAcqs.size(); i_acq++){

    #ifdef DEBUG 
    std::cout << " - Start with Acq #" << i_acq << std::endl;
    #endif

    TCT::acquisition_single* acq = &AllAcqs[i_acq];
    ana.AcqsAnalyser(acq, i_acq, &AcqAvg);

    #ifdef DEBUG 
    std::cout << *acq << std::endl;
    #endif

    if( ana.AcqsSelecter(acq) ) Nselected++;
    ana.AcqsProfileFiller(acq, &AcqAvg);

  }

  /*std::cout << AllAcqs.size() << std::endl;
    for(uint32_t i_acq = 0; i_acq < AllAcqs.size(); i_acq++){
    std::cout << " - Start with Acq #" << i_acq << std::endl;

    acq = AllAcqs[i_acq];
    std::cout << acq << std::endl;

    }*/
  ana.AcqsWriter(&dummyDUT2, &AllAcqs, &AcqAvg);

  std::cout << "Nselected = " << Nselected << std::endl;
  std::cout << "ratio of selected acqs = " << (float)Nselected/AllAcqs.size()*100. << "%" << std::endl;



  //theApp.Run(kTRUE); 
  //char key = getchar();

  std::cout << "end " << PACKAGE_NAME << std::endl;
  //theApp.Terminate();
  return 0;
}
