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

  TCT::util ana_card;

  std::string proj_folder = "def";
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-r")) {
      proj_folder = argv[++i];
      std::cout << "Root folder of TCT project is " << proj_folder << std::endl;
    }
    if (!strcmp(argv[i],"-af")) {
      std::ifstream ana_file (argv[++i]);
      std::cout << " Try to read analysis file from " << argv[i] << std::endl;
      ana_card.parse(ana_file);
    }
    // !! add check for certain vital options, if not passed, break
  }



  // create analysis object to stear the analysis. only one at the moment, no vector of objects. should be sufficient if analysis parameters are the same for every analysed folder
  TCT::analysis ana(ana_card.ID_val());
  //std::cout << ana << std::endl;

  for( auto i : ana_card.ID_val()){
    if(i.first == "ProjectFolder") {
      if(proj_folder != "def") std::cout << " Project folder from command line overwritten by value from analysis card. " << std:: endl;
      proj_folder = i.second;
    }
  }

  if(proj_folder == "def") {
    std::cout << "   *** Project folder not specified, neither in command line nor in analysis file!\n\n   ***STOPPING" << std::endl;
    exit(1);

  }

  if(!ana_card.IsRead()) {
    std::cout << "   *** No analysis card passed, please use -af option and specify file " << std::endl;
    exit(1);
  }

  // create sample object from sample card
  if(ana.SampleCard() == "def") {
    std::cout << "   *** No SampleCard was passed, check your analysis card, if \"SampleCard = ...\" is specified " << std::endl;
    exit(1);
  } else std::cout << "\n Try to read sample card from : " << ana.SampleCard() << std::endl;

  TCT::util sample_card;
  std::ifstream sample_file (ana.SampleCard());
  sample_card.parse(sample_file);

  if(sample_card.ID_val().size() == 0) {
    std::cout << "   *** No or only empty sample card found. Check that \"SampleCard = ...\" is specified correctly in analysis card.\n Then check identifiers and values in sample card" << std::endl;
    exit(1);
  }

  TCT::sample sample(sample_card.ID_val());
  //std::cout << sample << std::endl;	// print basic parameters of the sample

  if(sample.SampleID() == "def") {
    std::cout << "   *** No SampleID specified in sample card. Check that \"SampleID = ...\" is specified" << std::endl;
    exit(1);
  }

  // find all subfolders in datafolder
  if(ana.DataFolder() == "def") {
    std::cout << "   * No data folder was specified in analysis card. Check your analysis card, that \"DataFolder = ...\" is specified correctly" << std::endl;
    exit(1);
  }
  std::cout << "The sample's data folder is = " << ana.DataFolder() << ", searching data in subfolder(s) " <<  std::endl;
  std::map<std::string, std::vector<std::string>> folder_struc;
  std::vector<std::string> dirs;
  std::vector<std::string> dirs2;	// push folder for every subfolder (linearisation of folder matrix)
  std::vector<std::string> subdirs2;	// push folder for every subfolder (linearisation of folder matrix)
  std::vector<std::string> pathndirs;
  std::vector<int> Nsubdirs;

  void *dirp = gSystem->OpenDirectory(ana.DataFolder().c_str());
  if (!dirp) {
    std::cout << "   *** Data Folder not found. Check \"DataFolder\" in analysis card." << std::endl;
    exit(1);
  }
  char *direntry;
  uint32_t counterdir   = 0;
  uint32_t countersubdir= 0;
  while ((direntry = (char*)gSystem->GetDirEntry(dirp))) {
    if ( strstr(direntry,"..") || strstr(direntry,".") ) continue;
    counterdir++;
    //std::cout << "- " << direntry << std::endl;
    dirs.push_back((std::string)direntry);


    std::string subDataFolder	= ana.DataFolder() + "/" + direntry;
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
      std::string fullpath = ana.DataFolder() + "/" + direntry + "/" + subdirentry + "/";
      //std::cout << fullpath << std::endl;
      pathndirs.push_back(fullpath);

    }
    Nsubdirs.push_back(countersubdir);
    countersubdir = 0;

    folder_struc[dirs[counterdir-1]] = subdirs;


  }

  for (auto i : Nsubdirs) {
    countersubdir +=i;
    //std::cout << " i = " << i << std::endl;
  }

  std::cout << " Found the following folder structure: " << std::endl;
  for(auto i : folder_struc) {
    for(auto j : i.second)
      std::cout << i.first << " " <<  j << " " << "\n";
  }

  std::cout << " In total, found " << countersubdir << " subfolder(s) " << std::endl;





  //if(countersubdir > 0) {
  uint32_t counter = 0;
  while(1){


#ifdef DEBUG 
    std::cout << " Start with subfolder # " << counter << std::endl;
#endif
    //int i = 0; i < countersubdir; i++) {
    // create vec with acq_singles in it
    std::vector<TCT::acquisition_single> AllAcqs;

    if(countersubdir == 0) {
      // check if DataFolder() ends on "/", if not, add it
      std::string path = ana.DataFolder();
      if (path.length() > 0) {
	std::string::iterator it = path.end() - 1;
	if (*it != '/') {
	  path.append("/");
	}
      }
      ana.SetDataFolder(path);
      // push DataFolder to pathndirs (this is a hack to cope with Nsubdir == 0...)
      pathndirs.push_back(ana.DataFolder());
    }
    // create measurement object from one subdir for each cycle
    TCT::measurement meas(pathndirs[counter]);

    if(!meas.AcqsLoader(&AllAcqs, ana.MaxAcqs())) {
      std::cout << " Folder empty! Skipping folder" << std::endl; 
      continue;
    };

    // now create instance of avg acquisition using Nsamples from loaded files
    TCT::acquisition_avg AcqAvg(AllAcqs[0].Nsamples());
    AcqAvg.SetPolarity(AllAcqs[0].Polarity());

    //now analyse all acquisitions
    int Nselected = 0;

#ifdef DEBUG 
    std::cout << "Size of AllAcqs = " << AllAcqs.size() << std::endl;
#endif

    AcqAvg.SetNanalysed(AllAcqs.size());
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
      AcqAvg.SetNselected(Nselected);
      ana.AcqsProfileFiller(acq, &AcqAvg);

    }

    ana.SetOutSample_ID(sample.SampleID());
    if(countersubdir > 0) {
      ana.SetOutSubFolder(dirs2[counter]);
      ana.SetOutSubsubFolder(subdirs2[counter]);
    } else {
      ana.SetOutSubFolder("def");
      ana.SetOutSubsubFolder("def");
    }


    if(ana.SaveToFile()) 
      if(countersubdir > 0) ana.AcqsWriter(&AllAcqs, &AcqAvg);
      else ana.AcqsWriterNoSubs(&AllAcqs, &AcqAvg);

    std::cout << "   Nselected = " << Nselected << std::endl;
    std::cout << "   ratio of selected acqs = " << Nselected << " / " << AllAcqs.size() << " = " << (float)Nselected/AllAcqs.size()*100. << "%\n\n" << std::endl;

    // now take care of memory management
    // delete remaning TH1Fs in acquisition_single and then clear AllAcqs
    for(int j = 0; j < AllAcqs.size(); j++) {
      AllAcqs[j].Clear();
    }
    AllAcqs.clear();

    counter++;
    if( countersubdir == 0) break;
    if(counter == countersubdir) break;
#ifdef DEBUG 
    std::cout << "   AllAcqs has " << AllAcqs.size() << " objects left" << std::endl;
#endif

  } 

  //theApp.Run(kTRUE); 
  //char key = getchar();

  std::cout << "end " << PACKAGE_NAME << "\n" << std::endl;
  //theApp.Terminate();
  return 0;
  }
