/**
 * \file
 * \brief Main program code.
 */

//  includes from standard libraries
#include <iostream>
#include <regex>
#include <fstream>
#include <vector>
#include <map>

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
    std::cout	<< " No root folder specified. Please execute with path to projcet: > ./tct-analysis -r /home/<user>/<my-path>/TCT-analysis/" << std::endl;
    std::cout	<< " Or pass analysis file: > ./tct-analysis -af /home/<user>/<my-path>/TCT-analysis/testanalysis/ana.txt" << std::endl;
    std::cout	<< " Or like this: > ./tct-analysis -af ../testanalysis/ana.txt" << std::endl;
    std::cout	<< "\n Options are \n" 
      << "   -af <analysis file> (see sample analysis file for example)\n"
      << "   -r <project folder> (e.g. /home/<user>/<my-path>/TCT-analysis/ \n"
      << "   -sa (to save all single acquisition in root file (blows up root file)" // !! needs implementation
      << std::endl;
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

  // find all subfolders in datafolder
  std::string DataFolder	= proj_folder + "/testdata/S57";
  std::cout << "The sample folder is = " << DataFolder << ", searching data in subfolder(s) " <<  std::endl;
  std::map<std::string, std::vector<std::string>> folder_struc;
  std::vector<std::string> dirs;
  std::vector<std::string> dirs2;	// push folder for every subfolder (linearisation of folder matrix)
  std::vector<std::string> subdirs2;	// push folder for every subfolder (linearisation of folder matrix)
  std::vector<std::string> pathndirs;
  std::vector<int> Nsubdirs;

  void *dirp = gSystem->OpenDirectory(DataFolder.c_str());
  if (!dirp) {
    std::cout << "Data Folder not found" << std::endl;
    return 1;
  }
  char *direntry;
  uint32_t counterdir   = 0;
  uint32_t countersubdir= 0;
  while ((direntry = (char*)gSystem->GetDirEntry(dirp))) {
    if ( strstr(direntry,"..") || strstr(direntry,".") ) continue;
    counterdir++;
    //std::cout << "- " << direntry << std::endl;
    dirs.push_back((std::string)direntry);


    std::string subDataFolder	= DataFolder + "/" + direntry;
    //std::cout << "Subdir = " << subDataFolder << std::endl;
    std::vector<std::string> subdirs;
    void *subdirp = gSystem->OpenDirectory(subDataFolder.c_str());
    char *subdirentry;
    while ((subdirentry = (char*)gSystem->GetDirEntry(subdirp))) {
      if ( strstr(subdirentry,"..") || strstr(subdirentry,"." ) ) continue;
      countersubdir++;
      //std::cout << "-- " << subdirentry << std::endl;
      subdirs.push_back(subdirentry);
      dirs2.push_back((std::string)direntry);
      subdirs2.push_back((std::string)subdirentry);
      std::string fullpath = DataFolder + "/" + direntry + "/" + subdirentry + "/";
      //std::cout << fullpath << std::endl;
      pathndirs.push_back(fullpath);

    }
    Nsubdirs.push_back(countersubdir);
    countersubdir = 0;

    folder_struc[dirs[counterdir-1]] = subdirs;


  }
  for (auto i : Nsubdirs) {
    countersubdir +=i;
    std::cout << " i = " << i << std::endl;
  }

  std::cout << " Found the following folder structure: " << std::endl;
  for(auto i : folder_struc) {
    for(auto j : i.second)
      std::cout << i.first << " " <<  j << " " << "\n";
  }

  std::cout << "Found " << countersubdir << " subfolders " << std::endl;

  // create analysis object to stear the analysis. only one at the moment, no vector of objects. should be sufficient if analysis parameters are the same for every analysed folder
  TCT::analysis ana(util.ID_val());
  //std::cout << ana << std::endl;

  // Folder to read sensor card from. !! to be moved to analysis card
  std::string SensorFolder	= proj_folder + "/testSensor";

  // create smaple object from sample card
  TCT::sample dummyDUT2(SensorFolder);       // define DUT
  //std::cout << dummyDUT2 << std::endl;	// print basic parameters of the sample
  dummyDUT2.ReadSampleCard();	// read SampleCard and set parameters accordingly

  std::string sampleID = "S57"; // !! should come from sensor card
  dummyDUT2.SetSampleID(sampleID); // !! integrate into ReadSampleCard()
  //std::cout << dummyDUT2 << std::endl;


  for(int i = 0; i < countersubdir; i++) {
    // create vec with acq_singles in it
    std::vector<TCT::acquisition_single> AllAcqs;

    // create measurement object from one subdir for each cycle
    TCT::measurement meas(pathndirs[i]);

    if(!meas.AcqsLoader(&AllAcqs, MaxAcqs)) return 1;

    // now create instance of avg acquisition using Nsamples from loaded files
    TCT::acquisition_avg AcqAvg(AllAcqs[0].Nsamples());
    AcqAvg.SetPolarity(AllAcqs[0].Polarity());

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
      if(ana.DoSmearing()) ana.AcqsSmearer(acq, ana.AddNoise(), false);
      ana.AcqsAnalyser(acq, i_acq, &AcqAvg);
      if(ana.DoSmearing()) ana.AcqsSmearer(acq, false, ana.AddJitter()); // AcqsAnalyser removes jitter by determining each acqs delay. Hence, to add jitter, delay has to be manipulated after AcqsAnalyser (and before filling of profile

#ifdef DEBUG 
      std::cout << *acq << std::endl;
#endif

      if( ana.AcqsSelecter(acq) ) {
        Nselected++;
	acq->SetSelect(true);
      }
      ana.AcqsProfileFiller(acq, &AcqAvg);

    }

    ana.SetOutSample_ID(dummyDUT2.SampleID());
    ana.SetOutTemp(dirs2[i]);
    ana.SetOutVolt(subdirs2[i]);

    ana.AcqsWriter(&AllAcqs, &AcqAvg);

    std::cout << "   Nselected = " << Nselected << std::endl;
    std::cout << "   ratio of selected acqs = " << Nselected << " / " << AllAcqs.size() << " = " << (float)Nselected/AllAcqs.size()*100. << "%\n\n" << std::endl;

    // now take care of memory management
    // delete remaning TH1Fs in acquisition_single and then clear AllAcqs
    for(int j = 0; j < AllAcqs.size(); j++) {
    	AllAcqs[j].Clear();
    }
    AllAcqs.clear();

#ifdef DEBUG 
    std::cout << "   AllAcqs has " << AllAcqs.size() << " objects left" << std::endl;
#endif

  }

  //theApp.Run(kTRUE); 
  //char key = getchar();

  std::cout << "end " << PACKAGE_NAME << std::endl;
  //theApp.Terminate();
  return 0;
}
