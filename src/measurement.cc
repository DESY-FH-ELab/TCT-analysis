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

  bool measurement::AcqsLoader(std::vector<TCT::acquisition_single> *allAcqs, uint32_t MaxAcqs, bool LeCroyRAW){

    //if(debug) std::cout << "start PulseCheck" << std::endl;

    if(MaxAcqs == 0) {
      std::cout << "MaxAcqs is zero, stopping " << std::endl;
      return false;
    }
    if(MaxAcqs == -1) {
      std::cout << "MaxAcqs is -1, reading all files " << std::endl;
    }

    //std::cout << " MaxAcqs = " << MaxAcqs << std:: endl;

    const char* filedir = _DataInFolder.c_str(); // !! change to encapsulation 
    /*
     */
    //  gROOT->ProcessLine("#include <vector>");
    std::cout << " read files from: " << filedir << std::endl;
    //if (!AllTests) Init();
    //GetCuts(filedir);

    // get list of files in filedir
    void *dir = gSystem->OpenDirectory(filedir);
    const char *infile;
    uint32_t nfiles = 0;
    //std::cout << "do read in" << std::endl;
    //char key = getchar();
    std::cout<<"Parsing oscilloscope data"<<std::endl;
    if(LeCroyRAW) {
        std::cout << " Parsing data using LeCroy RAW reader " << std::endl;
        while((infile = gSystem->GetDirEntry(dir))) {
            if (strstr(infile,".trc")) {
                char pathandfile[250];
                strcpy(pathandfile,filedir);
                strcat(pathandfile,infile);
                if(nfiles < 3) std::cout << "  read file from: " << pathandfile << std::endl;
                if(nfiles == 3) std::cout << " suppressing further 'read from' info" << std::endl;

                TCT::acquisition_single acq(nfiles);
                std::string fullfname = pathandfile;
                bool read = acq.ReadRAW(fullfname, nfiles);
                nfiles++;
                //std::cout << "nfiles: " << nfiles;
                allAcqs->push_back(acq);
                if (nfiles > MaxAcqs-1) break;
            }
        }
    }
    else {
        std::cout << " Parsing data using *.txt reader " << std::endl;
        FILE *file;
        while((infile = gSystem->GetDirEntry(dir))) {
            if (strstr(infile,".txt") && !strstr(infile,".swp") ) {
                char pathandfile[250];
                strcpy(pathandfile,filedir);
                strcat(pathandfile,infile);
                if(nfiles < 3) std::cout << "  read file from: " << pathandfile << std::endl;
                if(nfiles == 3) std::cout << " suppressing further 'read from' info" << std::endl;
                file = fopen(pathandfile,"r");
                if (!file) {
                    std::cout << "   *** Can't open file! Exiting!" <<std::endl;
                    exit(1);
                }

                TCT::acquisition_single acq(nfiles);
                bool read = acq.Read(file, nfiles);
                fclose(file);
                nfiles++;
                //std::cout << "nfiles: " << nfiles;
                allAcqs->push_back(acq);
                if (nfiles > MaxAcqs-1) break;
            }
        }
    }

    gSystem->FreeDirectory(dir);



    if(nfiles == 0) { 
      std::cout << " -> no files read " << std::endl; return false;
    } else {
      std::cout << "   -> Found " << allAcqs->size() << " acquisitions, proceed with analysis" <<std::endl;
    }

    return true;
  } // end of LoadAcqs()


}
