/**
 * \file
 * \brief Implementation of TCT::Scanning methods.
 */

// STD includes
#include<string>
#include<sstream>

// TCT includes
#include "scanning.h"
#include "TCTReader.h"
#include "TCTModule.h"
#include "modules/ModuleLaserAnalysis.h"

// ROOT includes
#include "TH1F.h"
#include "TSystem.h"t

// External includes

namespace TCT {
/// Analysis manager
#ifndef USE_GUI
        bool Scanning::ReadTCT(char* filename, tct_config* config1) {
#else
/// Analysis manager
        bool Scanning::ReadTCT(char* filename, tct_config* config1, Ui::ConsoleOutput *progress) {
#endif
            /** Opens TCT data file, checks it, creates root file, runs analysis according to the config file.
             *  \param[in] filename Name of the TCT data file
             *  \param[in] config1 Pointer to the TCT::tct_config class
             *  \param[in] progress Pointer to the progress bar
              */
        config = config1;

        // -3 is the time shift, you can shift a signal to start at t=0. FIXME
        stct = new TCTReader(filename,-3,2);

        // Function corrects the baseline (DC offset) of all wafeforms
        // Float_t xc ; time denoting the start of the pulse
        //              correction factor is calculated from all the bins before xc
        // it integrates from first bin to bin with t=10, and then shifts by the mean value
        stct->CorrectBaseLine(config->CorrectBias());

        // CheckData: check if channels are set in config file 
        if(!CheckData()) {std::cout<<"File "<<filename<<" contains not enough data for selected operations. Skipping."<<std::endl; delete stct; return false;}

        //create output file
        CreateOutputFile();

        // write sample signals and separate waveforms
        Separate_and_Sample();
#ifdef USE_GUI
        if(config->FSeparateWaveforms()) progress->setValue(progress->value()+1);
#endif


        for(int i=0;i<config->GetNumberOfModules();i++) {
            if(config->GetModule(i)->isEnabled() && config->TCT_Mode()==(int)config->GetModule(i)->GetType()) {
                config->GetModule(i)->Do(stct,f_rootfile);
#ifdef USE_GUI
                progress->setValue(progress->value()+1);
#endif
            }
        }

        if(config->CH_PhDiode()) {
            ModuleLaserAnalysis* laser_analysis = new ModuleLaserAnalysis(config,"Laser_Analysis",_Top,"Analyse Laser Charge");
            laser_analysis->Do(stct,f_rootfile);
            delete laser_analysis;
        }

        f_rootfile->Close();

        delete stct;

        return true;
    }

/// Creates output file
bool Scanning::CreateOutputFile() {

    std::string outfolder	= config->OutFolder();
    std::string outpath  = outfolder + "/" + config->OutSample_ID();
    gSystem->MakeDirectory(outpath.c_str());

    char *InpName = stct->FileName;
    char * pch;
    pch = strtok (InpName,"/");
    while (pch != NULL){
        InpName = pch;
        pch = strtok (NULL, "/");
    }
    std::string inputname(InpName);
    std::string pathandfilename;
    pathandfilename = outpath  + "/" + config->OutSample_ID() + "_" + inputname + "_";
    char name[100];
    sprintf(name,"%02d.%02d.%02d",stct->Date[0],stct->Date[1],stct->Date[2]);
    pathandfilename = pathandfilename + name + ".root";

    std::cout << "Output file was created: " << pathandfilename << std::endl;
    f_rootfile = new TFile(pathandfilename.c_str(),"RECREATE","TCTanalyser");
    f_rootfile->cd();

    if(f_rootfile) return true;
    else {
        std::cout<<"!!! Failed to open root file!"<<std::endl;
        return false;
    }
}

/// Writes sample signals and separate waveforms
    bool Scanning::Separate_and_Sample() {

        TH1F *sample_hist;
        f_rootfile->mkdir("sample_signals");
        f_rootfile->cd("sample_signals");


        if(config->FSeparateWaveforms()) f_rootfile->mkdir("detector_signals");

        for(int i=0;i<4;i++) { // scan over channels
            if(stct->WFOnOff[i]) {
                sample_hist = stct->GetHA(i,0,0,0,0,0); // get one sample from each channel
                sample_hist->Write();
            }
            if(config->FSeparateWaveforms() && (i+1)==config->CH_Det() && stct->WFOnOff[i]) { //loop over all waveforms
                f_rootfile->cd("detector_signals");
                for(int l=0;l<stct->NU2;l++) {
                    for(int n=0;n<stct->NU1;n++) {
                        for(int j=0;j<stct->Nz;j++) {
                            for(int k=0;k<stct->Ny;k++) {
                                for(int m=0;m<stct->Nx;m++) {
                                    sample_hist = stct->GetHA(config->CH_Det()-1,m,k,j,n,l);
                                    sample_hist->Write();
                                }
                            }
                        }

                    }
                }
                f_rootfile->cd("sample_signals");
            }
        }
        f_rootfile->cd();

    }

/// Checks channels configuration
    bool Scanning::CheckData() {

        std::cout<<"Checking Channels set:"<<std::endl;
        std::cout<<"\t- Detector Signal Channel: "<< config->CH_Det() <<std::endl;
        if(config->CH_Det()) {
            if(stct->WFOnOff[config->CH_Det()-1]) std::cout<<"\t\t Data OK"<<std::endl;
            else {
                std::cout<<"\t\t This channel has no data. Please check the settings."<<std::endl;
                return false;
            }
        }
        else {
            std::cout<<"\t\t No data channel specified! Non-sense!"<<std::endl;
            return false;
        }
        std::cout<<"\t- Trigger Channel: "<< config->CH_Trig() <<std::endl;
        if(config->CH_Trig()) {
            if(stct->WFOnOff[config->CH_Trig()-1]) std::cout<<"\t\t Data OK"<<std::endl;
            else {
                std::cout<<"\t\t This channel has no data. Please check the settings."<<std::endl;
                return false;
            }
        }
        std::cout<<"\t- Photodiode Channel: "<< config->CH_PhDiode() <<std::endl;
        if(config->CH_PhDiode()) {
            if(stct->WFOnOff[config->CH_PhDiode()-1]) std::cout<<"\t\t Data OK"<<std::endl;
            else {
                std::cout<<"\t\t This channel has no data. Please check the settings."<<std::endl;
                return false;
            }
        }

        std::cout<<"Channel test passed. Processing..."<<std::endl;
        return true;
    }


}


