/**
 * \file
 * \brief Implementation of measurement methods
 */

#include<string>
#include <list>

#include "measurement.h"
//#include "acquisition.h"

//  includes from ROOT libraries
#include "TSystem.h"
#include "TFile.h"

namespace TCT {

  bool measurement::AcqsLoader(std::vector<TCT::acquisition_single> *allAcqs, uint32_t MaxAcqs){

    //if(debug) std::cout << "start PulseCheck" << std::endl;

    if(MaxAcqs == 0) {
      std::cout << "MaxAcqs is defaulted to zero, stopping " << std::endl;
      return false;
    }

    const char* filedir = _DataInFolder.c_str(); // !! change to encapsulation 
    /*
     */
    //  gROOT->ProcessLine("#include <vector>");
    printf("read files from %s\n",filedir);
    //if (!AllTests) Init();
    //GetCuts(filedir);

    // get list of files in filedir
    void *dir = gSystem->OpenDirectory(filedir);
    const char *infile;
    uint32_t nfiles = 0;
    //std::cout << "do read in" << std::endl;
    //char key = getchar();
    FILE *file;
    while(infile = gSystem->GetDirEntry(dir)) {
      if (strstr(infile,".txt") && !strstr(infile,".swp") ) {
	//printf("file: %s\n",infile);
	Char_t pathandfile[250];
	strcpy(pathandfile,filedir);
	strcat(pathandfile,infile);
	if(nfiles < 10) std::cout << "read file from: " << pathandfile << std::endl;
	if(nfiles == 10) std::cout << "suppressing further 'read from' info" << std::endl;
	file = fopen(pathandfile,"r");
	if (!file) {
	  std::cout <<"Can't open file"<<std::endl;
	  exit(1);
	}

	TCT::acquisition_single acq(nfiles); 
	//std::cout << acq->BiasVolt() << std::endl;
	//std::cout << "start Read() " << nfiles << std::endl;
	bool read = acq.Read(file, nfiles);
	fclose(file);
	//delete[] file; // doesnt work
	//if (read) std::cout << "File read!" << std::endl;
	nfiles++;
	//std::cout << "nfiles: " << nfiles;
	allAcqs->push_back(acq);
	if (nfiles > MaxAcqs && MaxAcqs > 0) break;
      }
    }
    gSystem->FreeDirectory(dir);
    //std::cout << "end read in" << std::endl;
    //char key = getchar();



    if(nfiles == 0) { 
      std::cout << " no files read, exiting..." << std::endl; return false;
    } else {
      std::cout << "Found " << allAcqs->size() << "acquision, proceed with analysis\n" <<std::endl;
    }

    return true;
  } // end of LoadAcqs()


}
