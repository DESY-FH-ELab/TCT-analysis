/**
 * \file
 * \brief Main program code.
 */

//  includes from standard libraries
#include <iostream>

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

  // !! for all file : 	AllTests->Add((BMTest *)pulse);	create vector of acqs

  std::string MeasFolder = "/home/hjansen/Diamond/data/S57/295k/500V/";
  std::string SensorFolder = "/home/hjansen/Sensors/testSensor/";

  TCT::sample dummyDUT;       // define DUT
  TCT::sample dummyDUT2(SensorFolder);       // define DUT
  std::cout << dummyDUT2 << std::endl;	// print basic parameters of the sample
  dummyDUT2.ReadSampleCard();	// read SampleCard and set parameters accordingly

  std::string sampleID = "DummySample1";
  dummyDUT2.SetSampleID(sampleID); 
  std::cout << dummyDUT2 << std::endl;

  TCT::acquisition_single acq0; 
  acq0.SetBiasVolt(600.);
  std::cout << acq0.BiasVolt() << std::endl;

  std::vector<TCT::acquisition_single> AllAcqs;
  //AllAcqs.push_back(acq0);
  
  //for (uint32_t i = 0; i < 4; i++{
  //  AllAcqs.push_back()
 // }
  // analyser
  //TCT::vdrift_ana vdrift;

  //TCT::param param;

  TCT::measurement meas(MeasFolder);
  std::cout << meas << std::endl;

  meas.AcqsLoader(AllAcqs);//, 4); // change to take parameter from param



  // print out method for analysers?!
  //std::cout << "dummy" << std::endl;


  //theApp.Run(kTRUE); 
  char key = getchar();

  std::cout << "end " << PACKAGE_NAME << std::endl;
  //theApp.Terminate();
  return 0;
}
