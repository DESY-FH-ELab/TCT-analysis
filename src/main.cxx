/**
 * \file
 * \brief Main program code.
 */

//  includes from standard libraries
#include <iostream>
#include <regex>

//  includes from TCT classes
#include "sample.h"
#include "config.h"
#include "param.h"
#include "acquisition.h"
#include "measurement.h"

//  includes from ROOT libraries
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

//#inculde "anaylyser.h" // inherits from class sample/measurement?

//using namespace TCT; // namespace of TCT_analysis is "TCT"
//using namespace std;

int main(int argc, char **argv)
{
  std::cout << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;

  //std::regex e("[0-9]");
  std::smatch sm; std::regex e ("(sub)(.*)");
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

  std::string DataFolder = "/home/hjansen/Diamond/data/S57/295k/500V/";
  std::string OutFolder = "/home/hjansen/Diamond/results";
  std::string SensorFolder = "/home/hjansen/Sensors/testSensor/";

  TCT::sample dummyDUT2(SensorFolder);       // define DUT
  std::cout << dummyDUT2 << std::endl;	// print basic parameters of the sample
  dummyDUT2.ReadSampleCard();	// read SampleCard and set parameters accordingly

  std::string sampleID = "S57";
  dummyDUT2.SetSampleID(sampleID); 
  std::cout << dummyDUT2 << std::endl;

  std::vector<TCT::acquisition_single> AllAcqs;


  //TCT::param param;

  TCT::measurement meas(DataFolder, OutFolder);
  std::cout << meas << std::endl;

  if(!meas.AcqsLoader(&AllAcqs)) return 1;//, 4); // change to take parameter from param
  // !! analyser functions belong to measurement class a.t.m., disentangle ?

  // now create instance of avg acquisition using Nsamples from loaded files
  TCT::acquisition_avg AcqAvg(AllAcqs[0].Nsamples());

  //now analysis all acquisitions
  int Nselected = 0;
  std::cout << AllAcqs.size() << std::endl;
  for(uint32_t i_acq = 0; i_acq < AllAcqs.size(); i_acq++){
    std::cout << " - Start with Acq #" << i_acq << std::endl;

    TCT::acquisition_single acq = AllAcqs[i_acq];
    std::cout << acq << std::endl;

    meas.AcqsAnalyser(&acq, i_acq, &AcqAvg);
    if( meas.AcqsSelecter(&acq) ) Nselected++;
    meas.AcqsProfileFiller(&acq, &AcqAvg);

  }

  meas.AcqsWriter(&dummyDUT2, &AllAcqs, &AcqAvg);

  std::cout << "Nselected = " << Nselected << std::endl;
  std::cout << "ratio of selected acqs = " << (float)Nselected/AllAcqs.size()*100. << "%" << std::endl;



  //theApp.Run(kTRUE); 
  //char key = getchar();

  std::cout << "end " << PACKAGE_NAME << std::endl;
  //theApp.Terminate();
  return 0;
}
