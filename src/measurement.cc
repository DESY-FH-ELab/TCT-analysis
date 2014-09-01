/**
 * \file
 * \brief Implementation of measurement methods
 */

#include<string>
#include "measurement.h"
#include "acquisition.h"

//  includes from ROOT libraries
#include "TSystem.h"

namespace TCT {

  void measurement::AcqsLoader(std::vector<TCT::acquisition_single> AllAcqs, uint32_t MaxAcqs){

    //if(debug) std::cout << "start PulseCheck" << std::endl;

    const char* filedir = _Folder.c_str();  
    /*
     */
    //  gROOT->ProcessLine("#include <vector>");
    printf("read files from %s\n",filedir);
    //if (!AllTests) Init();
    //GetCuts(filedir);

    // get list of files in filedir
    void *dir = gSystem->OpenDirectory(filedir);
    const char *infile;
    Int_t nfiles = 0;
    //std::cout << "do read in" << std::endl;
    //char key = getchar();
    FILE *file;
    while(infile = gSystem->GetDirEntry(dir)) {
      if (strstr(infile,".txt") && !strstr(infile,".swp") ) {
	//printf("file: %s\n",infile);
	Char_t pathandfile[250];
	strcpy(pathandfile,filedir);
	strcat(pathandfile,infile);
	if(nfiles < 5) std::cout << "read file from: " << pathandfile << std::endl;
	if(nfiles == 5) std::cout << "suppressing further 'read from' info" << std::endl;
	file = fopen(pathandfile,"r");
	if (!file) {
	  std::cout <<"Can't open file"<<std::endl;
	  exit(1);
	}

	TCT::acquisition_single acq; 
	//std::cout << acq.BiasVolt() << std::endl;
	//std::cout << "start Read() " << nfiles << std::endl;
	bool read = acq.Read(file, nfiles);
	fclose(file);
	//delete[] file; // doesnt work
	//if (read) std::cout << "File read!" << std::endl;
	nfiles++;
	//std::cout << "nfiles: " << nfiles;
	AllAcqs.push_back(acq);
	if (nfiles > MaxAcqs && MaxAcqs > 0) break;
      }
    }
    //std::cout << "end read in" << std::endl;
    //char key = getchar();

    std::cout << "Found " << AllAcqs.size() << "acquision, proceed with analysis\n" <<std::endl;

    
    if(nfiles == 0) { std::cout << " no files read, exiting..." << std::endl; return;}

    //std::cout << "delay cut: " << delay_cut << " delay win: " << delay_win << std::endl;
    //std::cout << " ampl cut: " << ampl_cut <<  "  ampl win: " << ampl_win << std::endl;
    //std::cout << "width cut: " << width_cut << " width win: " << width_win << std::endl;
    //if(cutRising) std::cout << "rise time cut: " << rise_cut << std::endl;

    //now analysis all pulses
    for(uint32_t i_acq = 0; i_acq < AllAcqs.size(); i_acq++){
      std::cout << " - Start with Acq #" << i_acq << std::endl;
      //(AllAcqs(i_acq)).BasicAnalysis();
    }

      return;
  } // end of LoadAcqs()

  //void measurement::AcqsAnalyser(std::vector<TCT::acquisition_single> AllAcqs){

  //}

  //void measurement::AcqsWriter(std::vector<TCT::acquisition_single> AllAcqs){

    /*

    // write results
    string outhelp;

    string outfolder = "/home/hjansen/Diamond/results"; // all

    //std::cout << "size " <<   filedirstring.size() << std::endl;
    string str1 = filedirstring.substr(0,filedirstring.size()-1);
    string outvolt = str1.substr(str1.rfind("/")+1);
    std::cout << outvolt << std::endl;

    string str2 = str1.substr(0,str1.rfind("/"));
    string outtemp = str2.substr(str2.rfind("/")+1);
    std::cout << outtemp << std::endl;

    string sample = "S57_23";

    outhelp = outfolder + "/" + sample + "/" + outtemp + "/" + sample + "_" + outtemp + "_" + outvolt + ".root";
    std::cout << "outfile name: " << outhelp << std::endl;


    //WriteRoot(outhelp.c_str());
    if(nvolt > 1 || ntemp > 1){
      std::cout << "		clearing root" << std::endl;
      ClearRoot();
    }


    if(debug) std::cout << "end PulseCheck" << std::endl;

    */
    //}

}
