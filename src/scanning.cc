/**
 * \file
 * \brief Implementation of daq reading methods
 */

// STD includes
#include<string>
#include<sstream>

// TCT includes
#include "scanning.h"

// ROOT includes
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "TImage.h"
#include "TDirectory.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TLegend.h"
// External includes
//#include "TCTScan.h"
#include "TCTReader.h"

namespace TCT {
    bool Scanning::ReadTCT(char *filename, tct_config *config1) {
        config = config1;

        // -3 is the time shift, you can shift a signal to start at t=0. FIXME
        stct = new TCTReader(filename,-3,2);

        // Function corrects the baseline (DC offset) of all wafeforms
        // Float_t xc ; time denoting the start of the pulse
        //              correction factor is calculated from all the bins before xc
        // it integrates from first bin to bin with t=10, and then shifts by the mean value
        stct->CorrectBaseLine(10.);

        // CheckData: check if channels are set in config file 
        if(!CheckData()) {std::cout<<"File "<<filename<<" contains not enough data for selected operations. Skipping."<<std::endl; return false;}

        std::string outfolder	= config->OutFolder();
        std::string outpath  = outfolder + "/" + config->OutSample_ID();
        gSystem->MakeDirectory(outpath.c_str());

        std::string pathandfilename;
        pathandfilename = outpath  + "/" + config->OutSample_ID() + "_";
        char name[100];
        sprintf(name,"%02d.%02d.%02d_%02d.%02d.%02d",stct->Date[0],stct->Date[1],stct->Date[2],stct->Date[3],stct->Date[4],stct->Date[5]);
        pathandfilename = pathandfilename + name + ".root";

        std::cout << "Output file was created: " << pathandfilename << std::endl;
        f_rootfile = new TFile(pathandfilename.c_str(),"RECREATE","TCTanalyser");
        f_rootfile->cd();

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
        if(config->DO_focus() && config->TCT_Mode()==0) DoTopFocus();
        if(config->DO_focus() && config->TCT_Mode()==1) DoEdgeFocus();
        //if(config->DO_focus() && config->TCT_Mode()==2) BottomDoFocus(); // FIXME needs implementation
        if(config->DO_EdgeDepletion() && config->TCT_Mode()==1) DoEdgeDepletion();
        if(config->DO_EdgeVelocity() && config->TCT_Mode()==1) DoEdgeVelocity();
        if(config->CH_PhDiode()) LaserPowerDrop();
        if(config->CH_PhDiode()) BeamSigma();

        f_rootfile->Close();

        delete stct;

        return true;
    }

    bool Scanning::DoTopFocus() {

        if(!CheckFocus()) {std::cout<<"No data for focusing. Skipping..."<<std::endl; return false;}
        TDirectory *dir_fsearch = f_rootfile->mkdir("FocusSearch");
        TDirectory *dir_fsearch_normed;
        if(config->CH_PhDiode()) dir_fsearch_normed = f_rootfile->mkdir("FocusSearch_Normed");
        dir_fsearch->cd();
        Int_t numO,numS;

        Float_t Ss,Os;                       // Optical axis step
        Float_t Sc0,Opt0;

        Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
        Int_t optic_axis=(config->OptAxis())-1;            // select optic axis (0=x,1=y,2=z)
        Int_t scanning_axis=config->ScAxis()-1;         // select scanning axis (0=x,1=y,2=z)

        SwitchAxis(optic_axis,numO,Os,Opt0);

        TGraph **cc = new TGraph*[numO];                  // charge collection graph
        TGraph **cc_norm;                  // charge collection graph
        TGraph **ph_charge = new TGraph*[numO];

        Float_t *width 		= new Float_t[numO];
        Float_t *pos 		= new Float_t[numO];
        Float_t *abs_pos 	= new Float_t[numO];
        Float_t *width_normed 	= new Float_t[numO];
        Float_t *pos_normed 	= new Float_t[numO];
        Float_t *strip_w 	= new Float_t[numO];
        Float_t *strip_w_normed = new Float_t[numO];
        Float_t *minQ 		= new Float_t[numO];
        Float_t *minQ_normed 	= new Float_t[numO];
        Float_t *optical_axis_co= new Float_t[numO];

        SwitchAxis(scanning_axis,numS,Ss,Sc0);

        //calculate the arb charge from the current
        CalculateCharges(ChNumber,optic_axis,numO,scanning_axis,numS,cc,config->FTlow(),config->FThigh());

        //laser charge distribution
        if(config->CH_PhDiode()) {

            //calculate the photodetector charge
            CalculateCharges(config->CH_PhDiode()-1,optic_axis,numO,scanning_axis,numS,ph_charge,config->FDLow(),config->FDHigh());

            dir_fsearch_normed->cd();
            ChargeCorrelationHist(cc,ph_charge,numO);
            dir_fsearch->cd();
        }

        // fit definition

        Float_t FWHM=config->FFWHM();               //   expected FWHM, can be specified in config file, dafaulted to 10. in header

        TF1 *ff2=new TF1("ff2","[2]/2.*(TMath::Erfc((x-[0])/[1]) + (TMath::Erf((x-[0]-[3])/[1]) + 1))",0,Ss*numS);

        double *yy_temp;
        Int_t i_max, i_min;
        Float_t max, min;

        for(int j=0;j<numO;j++){

            yy_temp = cc[j]->GetY();
            max = yy_temp[0];
            min = yy_temp[0];

            for(int i=1;i<numS;i++) {
                if(yy_temp[i]>max) {
                    max=yy_temp[i];
                    i_max=i;
                }
                if(yy_temp[i]<min) {
                    min=yy_temp[i];
                    i_min=i;
                }
            }

            SetFitParameters(ff2,i_min*Ss,FWHM,max,FWHM);

            cc[j]->Fit("ff2","Rq");
            gStyle->SetOptFit(1);
            width[j]=ff2->GetParameter(1)*2.35/TMath::Sqrt(2);
            pos[j]=ff2->GetParameter(0)+ ff2->GetParameter(3)/2;
            minQ[j] = ff2->Eval(ff2->GetParameter(0) + ff2->GetParameter(3)/2)/ff2->GetParameter(2);
            strip_w[j] = ff2->GetParameter(3);
            optical_axis_co[j]=Opt0+j*Os;
        }


        //calculating the normed charge distribution
        if(config->CH_PhDiode()) {

            cc_norm = NormedCharge(cc,ph_charge,numO);

            for(int j=0;j<numO;j++) {
                yy_temp = cc_norm[j]->GetY();
                max = yy_temp[0];
                min = yy_temp[0];

                for(int i=1;i<numS;i++) {
                    if(yy_temp[i]>max) {
                        max=yy_temp[i];
                        i_max=i;
                    }
                    if(yy_temp[i]<min) {
                        min=yy_temp[i];
                        i_min=i;
                    }
                }
                SetFitParameters(ff2,i_min*Ss,FWHM,max,FWHM);

                cc_norm[j]->Fit("ff2","Rq");
                gStyle->SetOptFit(1);

                width_normed[j]=ff2->GetParameter(1)*2.35/TMath::Sqrt(2);
                pos_normed[j]=ff2->GetParameter(0)+ ff2->GetParameter(3)/2;
                minQ_normed[j] = ff2->Eval(ff2->GetParameter(0) + ff2->GetParameter(3)/2)/ff2->GetParameter(2);
                strip_w_normed[j] = ff2->GetParameter(3);
            }
        }

        if(config->FSeparateCharges()) {
            dir_fsearch->cd();
            GraphSeparate(numO,cc,"charges","scanning distance, [#mum]","Charge,[arb.]","Charge vs Distance","Z = ",optical_axis_co);
            if(config->CH_PhDiode()) {
                dir_fsearch_normed->cd();
                GraphSeparate(numO,ph_charge,"photodiode","scanning distance, [#mum]","Charge,[arb.]","Laser Charge vs Distance","Z = ",optical_axis_co);
                GraphSeparate(numO,cc_norm,"charges_normed","scanning distance, [#mum]","Charge,[arb.]","Charge vs Distance","Z = ",optical_axis_co);
                dir_fsearch->cd();
            }
        }

        //plot graphs at different positions along optical axis
        MultiGraphWriter(numO,cc,"scanning distance [#mum]","Charge [arb.]","Charge Vs Distance","ChargeVsDistance");

        //plot graphs at different positions along optical axis with normed data
        if(config->CH_PhDiode()) {
            dir_fsearch_normed->cd();
            MultiGraphWriter(numO,cc_norm,"scanning distance [#mum]","Charge [arb.]","Charge Vs Distance","ChargeVsDistance_Normed");
            dir_fsearch->cd();
        }

        TF1 *ff_pol1 = new TF1("my_pol1","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol1->SetParName(0,"Focus Position");
        TF1 *ff_pol2 = new TF1("my_pol2","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol2->SetParName(0,"Focus Position");
        TF1 *ff_pol3 = new TF1("my_pol3","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol3->SetParName(0,"Focus Position");
        TF1 *ff_pol4 = new TF1("my_pol4","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol4->SetParName(0,"Focus Position");
        TPaveStats *st;

        //draw the gaussian beam profile
        TGraph *FWHMg = GraphBuilder(numO,optical_axis_co,width,"optical distance [#mum]","FWHM [#mum]","Gaussian Beam Profile");
        ff_pol1->SetParameter(0,FWHMg->GetMean());
        ff_pol1->SetParameter(2,FWHMg->GetMinimum());
        FWHMg->Fit("my_pol1","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg->Write("FWHM");

        // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance

        TGraph *MINQg=GraphBuilder(numO,optical_axis_co,minQ,"optical distance [#mum]","min Q[rel. to max]","Minimum Charge");
        ff_pol3->SetParameter(0,MINQg->GetMean());
        ff_pol3->SetParameter(2,MINQg->GetMinimum());
        MINQg->Fit("my_pol3","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)MINQg->FindObject("stats");
        st->SetFitFormat(".5g");
        MINQg->Write("MinCharge");

        // Find the missalignment between z and optical axis
        for(int j=0;j<numO;j++) abs_pos[j]=pos[j]+Sc0;

        TGraph *POSg=GraphBuilder(numO,optical_axis_co,abs_pos,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
        POSg->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg->Write("Missalignment");

        // Plotting best strip width
        GraphBuilder(numO,optical_axis_co,strip_w,"optical distance [#mum]","strip width [#mum]","Strip Width","StripWidth");

        if(config->CH_PhDiode()) {
            dir_fsearch_normed->cd();

            //draw the gaussian beam profile with normed charge
            TGraph *FWHMg1 = GraphBuilder(numO,optical_axis_co,width_normed,"optical distance [#mum]","FWHM [#mum]","Gaussian Beam Profile");
            ff_pol2->SetParameter(0,FWHMg1->GetMean());
            ff_pol2->SetParameter(2,FWHMg1->GetMinimum());
            FWHMg1->Fit("my_pol2","q");
            gStyle->SetOptFit(1);
            st = (TPaveStats*)FWHMg1->FindObject("stats");
            st->SetFitFormat(".5g");
            FWHMg1->Write("FWHM_Normed");

            // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance with normed data

            TGraph *MINQg1=GraphBuilder(numO,optical_axis_co,minQ_normed,"optical distance [#mum]","min Q[rel. to max]","Minimum Charge");
            ff_pol4->SetParameter(0,MINQg1->GetMean());
            ff_pol4->SetParameter(2,MINQg1->GetMinimum());
            MINQg1->Fit("my_pol4","q");
            gStyle->SetOptFit(1);
            st = (TPaveStats*)MINQg->FindObject("stats");
            st->SetFitFormat(".5g");
            MINQg1->Write("MinCharge_Normed");

            // Find the missalignment between z and optical axis
            for(int j=0;j<numO;j++) abs_pos[j]=pos_normed[j]+Sc0;

            TGraph *POSg1=GraphBuilder(numO,optical_axis_co,abs_pos,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
            POSg1->Fit("pol1","q");
            gStyle->SetOptFit(1);
            POSg1->Write("Missalignment_Normed");

            // Plotting best strip width
            GraphBuilder(numO,optical_axis_co,strip_w,"optical distance [#mum]","strip width [#mum]","Strip Width","StripWidth_Normed");

            dir_fsearch->cd();

        }

        delete width;
        delete pos;
        delete abs_pos;
        delete width_normed;
        delete pos_normed;
        delete strip_w;
        delete strip_w_normed;
        delete minQ;
        delete minQ_normed;
        delete optical_axis_co;

        return true;
    }

    bool Scanning::DoEdgeFocus() {

        if(!CheckFocus()) {std::cout<<"No data for focusing. Skipping..."<<std::endl; return false;}
        TDirectory *dir_fsearch = f_rootfile->mkdir("FocusSearch");
        TDirectory *dir_fsearch_normed;
        if(config->CH_PhDiode()) dir_fsearch_normed = f_rootfile->mkdir("FocusSearch_Normed");
        dir_fsearch->cd();

        Int_t numO,numS;

        Float_t Ss,Os;                       // Optical axis step
        Float_t Sc0,Opt0;

        Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
        Int_t optic_axis=(config->OptAxis())-1;            // select optic axis (0=x,1=y,2=z)
        Int_t scanning_axis=config->ScAxis()-1;         // select scanning axis (0=x,1=y,2=z)

        SwitchAxis(optic_axis,numO,Os,Opt0);

        TGraph **cc = new TGraph*[numO];                  // charge collection graph
        TGraph **cc_norm;                  // charge collection graph
        TGraph **ph_charge = new TGraph*[numO];

        Float_t *width_left = new Float_t[numO];
        Float_t *pos_left = new Float_t[numO];
        Float_t *abs_pos_left = new Float_t[numO];
        Float_t *width_right = new Float_t[numO];
        Float_t *pos_right = new Float_t[numO];
        Float_t *abs_pos_right = new Float_t[numO];
        Float_t *width_left_normed = new Float_t[numO];
        Float_t *pos_left_normed = new Float_t[numO];
        Float_t *width_right_normed = new Float_t[numO];
        Float_t *pos_right_normed = new Float_t[numO];
        Float_t *sensor_thick = new Float_t[numO];
        Float_t *sensor_thick_normed = new Float_t[numO];
        Float_t *optical_axis_co = new Float_t[numO];

        // set the scanning axis
        SwitchAxis(scanning_axis,numS,Ss,Sc0);

        //calculate the arb charge from the current
        CalculateCharges(ChNumber,optic_axis,numO,scanning_axis,numS,cc,config->FTlow(),config->FThigh());

        //laser charge distribution
        if(config->CH_PhDiode()) {

            //calculate the photodetector charge
            CalculateCharges(config->CH_PhDiode()-1,optic_axis,numO,scanning_axis,numS,ph_charge,config->FDLow(),config->FDHigh());

            dir_fsearch_normed->cd();
            ChargeCorrelationHist(cc,ph_charge,numO);
            dir_fsearch->cd();
        }

        //fit edges
        FindEdges(cc,numO,numS,Ss,pos_left,width_left,pos_right,width_right);
        for(int j=0;j<numO;j++) {
            sensor_thick[j]=pos_right[j]-pos_left[j];
            optical_axis_co[j]=Opt0+j*Os;
        }

        //calculating the normed charge distribution
        if(config->CH_PhDiode()) {

            cc_norm = NormedCharge(cc,ph_charge,numO);
            //find edges
            FindEdges(cc_norm,numO,numS,Ss,pos_left_normed,width_left_normed,pos_right_normed,width_right_normed);
            for(int j=0;j<numO;j++) sensor_thick_normed[j]=pos_right_normed[j]-pos_left_normed[j];
        }

        if(config->FSeparateCharges()) {
            dir_fsearch->cd();
            GraphSeparate(numO,cc,"charges","scanning distance, [#mum]","Charge,[arb.]","Charge vs Distance","Z = ",optical_axis_co);
            if(config->CH_PhDiode()) {
                dir_fsearch_normed->cd();
                GraphSeparate(numO,ph_charge,"photodiode","scanning distance, [#mum]","Charge,[arb.]","Laser Charge vs Distance","Z = ",optical_axis_co);
                GraphSeparate(numO,cc_norm,"charges_normed","scanning distance, [#mum]","Charge,[arb.]","Charge vs Distance","Z = ",optical_axis_co);
                dir_fsearch->cd();
            }
        }

        //plot graphs at different positions along optical axis
        MultiGraphWriter(numO,cc,"scanning distance [#mum]","Charge [arb.]","Charge Vs Distance","ChargeVsDistance");

        //plot graphs at different positions along optical axis with normed data
        if(config->CH_PhDiode()) {
            dir_fsearch_normed->cd();
            MultiGraphWriter(numO,cc_norm,"scanning distance [#mum]","Charge [arb.]","Charge Vs Distance","ChargeVsDistance_Normed");
            dir_fsearch->cd();
        }

        TF1 *ff_pol1 = new TF1("my_pol1","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol1->SetParName(0,"Focus Position");
        TF1 *ff_pol2 = new TF1("my_pol2","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol2->SetParName(0,"Focus Position");
        TPaveStats *st;

        TGraph *FWHMg_Left = GraphBuilder(numO,optical_axis_co,width_left,"Voltage [V]","FWHM [#mum]","Gaussian Beam Profile");
        ff_pol1->SetParameter(0,FWHMg_Left->GetMean());
        ff_pol1->SetParameter(2,FWHMg_Left->GetMinimum());
        FWHMg_Left->Fit("my_pol1","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg_Left->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg_Left->Write("FWHM_Left");

        TGraph *FWHMg_Right = GraphBuilder(numO,optical_axis_co,width_right,"Voltage [V]","FWHM [#mum]","Gaussian Beam Profile");
        ff_pol1->SetParameter(0,FWHMg_Right->GetMean());
        ff_pol1->SetParameter(2,FWHMg_Right->GetMinimum());
        FWHMg_Right->Fit("my_pol1","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg_Right->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg_Right->Write("FWHM_Right");

        // Find the missalignment between z and optical axis
        for(int j=0;j<numO;j++) {
            abs_pos_left[j]=pos_left[j]+Sc0;
            abs_pos_right[j]=pos_right[j]+Sc0;
        }

        TGraph *POSg_Left = GraphBuilder(numO,optical_axis_co,abs_pos_left,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
        POSg_Left->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg_Left->Write("Missalignment_Left");

        TGraph *POSg_Right = GraphBuilder(numO,optical_axis_co,abs_pos_right,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
        POSg_Right->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg_Right->Write("Missalignment_Right");

        // Plotting the sensor thickness
        GraphBuilder(numO,optical_axis_co,sensor_thick,"optical distance [#mum]","sensor thickness [#mum]","Sensor Thickness","SensorThickness");


        if(config->CH_PhDiode()) {
            dir_fsearch_normed->cd();

            //draw the gaussian beam profile with normed charge
            TGraph *FWHMg1_Left = GraphBuilder(numO,optical_axis_co,width_left_normed,"Voltage [V]","FWHM [#mum]","Gaussian Beam Profile");
            ff_pol2->SetParameter(0,FWHMg1_Left->GetMean());
            ff_pol2->SetParameter(2,FWHMg1_Left->GetMinimum());
            FWHMg1_Left->Fit("my_pol2","q");
            gStyle->SetOptFit(1);
            st = (TPaveStats*)FWHMg1_Left->FindObject("stats");
            st->SetFitFormat(".5g");
            FWHMg1_Left->Write("FWHM_Left_Normed");

            //draw the gaussian beam profile with normed charge
            TGraph *FWHMg1_Right = GraphBuilder(numO,optical_axis_co,width_right_normed,"Voltage [V]","FWHM [#mum]","Gaussian Beam Profile");
            ff_pol2->SetParameter(0,FWHMg1_Right->GetMean());
            ff_pol2->SetParameter(2,FWHMg1_Right->GetMinimum());
            FWHMg1_Right->Fit("my_pol2","q");
            gStyle->SetOptFit(1);
            st = (TPaveStats*)FWHMg1_Right->FindObject("stats");
            st->SetFitFormat(".5g");
            FWHMg1_Right->Write("FWHM_Right_Normed");

            // Find the missalignment between z and optical axis
            for(int j=0;j<numO;j++) {
                abs_pos_left[j]=pos_left_normed[j]+Sc0;
                abs_pos_right[j]=pos_right_normed[j]+Sc0;
            }

            TGraph *POSg1_Left = GraphBuilder(numO,optical_axis_co,abs_pos_left,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
            POSg1_Left->Fit("pol1","q");
            gStyle->SetOptFit(1);
            POSg1_Left->Write("Missalignment_Left_Normed");

            TGraph *POSg1_Right = GraphBuilder(numO,optical_axis_co,abs_pos_right,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
            POSg1_Right->Fit("pol1","q");
            gStyle->SetOptFit(1);
            POSg1_Right->Write("Missalignment_Right_Normed");

            // Plotting the sensor thickness
            GraphBuilder(numO,optical_axis_co,sensor_thick_normed,"optical distance [#mum]","sensor thickness [#mum]","Sensor Thickness","SensorThickness_Normed");
            dir_fsearch->cd();
        }

        delete width_left;
        delete pos_left;
        delete abs_pos_left;
        delete width_right;
        delete pos_right;
        delete abs_pos_right;
        delete width_left_normed;
        delete pos_left_normed;
        delete width_right_normed;
        delete pos_right_normed;
        delete sensor_thick;
        delete sensor_thick_normed;
        delete optical_axis_co;

        return true;
    }

    bool Scanning::DoEdgeDepletion() {

        if(!CheckEdgeDepletion()) {std::cout<<"No data for depletion voltage search. Skipping..."<<std::endl; return false;}
        TDirectory *dir_depl = f_rootfile->mkdir("DepletionVoltage");
        TDirectory *dir_depl_normed;
        if(config->CH_PhDiode()) dir_depl_normed = f_rootfile->mkdir("DepletionVoltage_Normed");
        dir_depl->cd();
        Int_t numVolt,numS;
        Float_t Ss,Sc0;

        Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
        Int_t volt_source=(config->VoltSource());            // select
        Int_t scanning_axis=config->ScAxis()-1;         // select scanning axis (0=x,1=y,2=z)
        Float_t *voltages;

        switch(volt_source)
          {
          case 1: numVolt=stct->NU1; voltages=stct->U1.fArray; break; //u1
          case 2: numVolt=stct->NU2; voltages=stct->U2.fArray; break; //u2
          }

        TGraph **cc = new TGraph*[numVolt];                  // charge collection graph
        TGraph **cc_norm = new TGraph*[numVolt];                  // charge collection graph
        TGraph **ph_charge = new TGraph*[numVolt];

        SwitchAxis(scanning_axis,numS,Ss,Sc0);

        //calculate charge profiles for different voltages
        CalculateCharges(ChNumber,volt_source+2,numVolt,scanning_axis,numS,cc,config->FTlow(),config->FThigh());

        //calculate laser charge profiles
        if(config->CH_PhDiode()) CalculateCharges(config->CH_PhDiode()-1,volt_source+2,numVolt,scanning_axis,numS,ph_charge,config->FDLow(),config->FDHigh());

        //calculating the normed charge distribution
        if(config->CH_PhDiode()) cc_norm=NormedCharge(cc,ph_charge,numVolt);

        // integrate charge through the detector

        Double_t right_edge,left_edge;

        FindEdges(cc[numVolt-1],numS,Ss,left_edge,right_edge);
        if(abs((right_edge-left_edge)-config->SampleThickness())>0.1*config->SampleThickness()) {
            std::cout<<"\tSorry, the detector thickness found for the highest bias voltage is more than 10% differes from the given in the configuration file. Most possible, that the detector is misaligned."<<std::endl;
            std::cout<<"\tLeft Edge = "<<left_edge<<" Right Edge = "<<right_edge<<std::endl;
            std::cout<<"\tOriginal Thickness is "<<config->SampleThickness()<<" micrometers"<<std::endl;
            std::cout<<std::endl;
        }
        else {
            std::cout<<"\tFound detector thickness is: "<<right_edge-left_edge<<" micrometers"<<std::endl;
            std::cout<<"\tIntegrating in the range "<<left_edge<<" - "<<right_edge<<std::endl;
            std::cout<<std::endl;
        }

        Float_t* total_charge = new Float_t[numVolt];
        //std::cout<<"charge left: "<<charge_left<<" charge right: "<<charge_right<<std::endl;
        for(int j=0;j<numVolt;j++) {
            total_charge[j] = GraphIntegral(cc[j],left_edge,right_edge);
        }
        if(total_charge[numVolt-1]<0) {
            for(int j=0;j<numVolt;j++) total_charge[j]=-total_charge[j];
        }

        FindEdges(cc_norm[numVolt-1],numS,Ss,left_edge,right_edge);
        if(abs((right_edge-left_edge)-config->SampleThickness())>0.1*config->SampleThickness()) {
            std::cout<<"\tNORMED Sorry, the detector thickness found for the highest bias voltage is more than 10% differes from the given in the configuration file. Most possible, that the detector is misaligned."<<std::endl;
            std::cout<<"\tNORMED Left Edge = "<<left_edge<<" Right Edge = "<<right_edge<<std::endl;
            std::cout<<"\tNORMED Original Thickness is "<<config->SampleThickness()<<" micrometers"<<std::endl;
            std::cout<<"\tNORMED Aborting the depletion voltage search"<<std::endl;
            return false;
        }
        else {
            std::cout<<"\tNORMED Found detector thickness is: "<<right_edge-left_edge<<" micrometers"<<std::endl;
            std::cout<<"\tNORMED Integrating in the range "<<left_edge<<" - "<<right_edge<<std::endl;
        }

        Float_t* total_charge_normed = new Float_t[numVolt];
        for(int j=0;j<numVolt;j++) {
            total_charge_normed[j] = GraphIntegral(cc_norm[j],left_edge,right_edge);
        }
        if(total_charge_normed[numVolt-1]<0) {
            for(int j=0;j<numVolt;j++) total_charge_normed[j]=-total_charge_normed[j];
        }

        //write separate charges
        if(config->FSeparateCharges()) {
            dir_depl->cd();
            GraphSeparate(numVolt,cc,"charges","scanning distance, [#mum]","Charge,[arb.]","Charge vs Distance","U = ",voltages);
            if(config->CH_PhDiode()) {
                dir_depl_normed->cd();
                GraphSeparate(numVolt,ph_charge,"photodiode","scanning distance, [#mum]","Charge,[arb.]","Laser Charge vs Distance","U = ",voltages);
                GraphSeparate(numVolt,cc_norm,"charges_normed","scanning distance, [#mum]","Charge,[arb.]","Charge vs Distance","U = ",voltages);
                dir_depl->cd();
            }
        }

        //plot graphs for different voltages
        MultiGraphWriter(numVolt,cc,"scanning distance [#mum]","Charge [arb.]","Charge Vs Distance","ChargeVsDistance");

        //plot graphs for different voltages with normed data
        if(config->CH_PhDiode()) {
            dir_depl_normed->cd();
            MultiGraphWriter(numVolt,cc_norm,"scanning distance [#mum]","Charge [arb.]","Charge Vs Distance","ChargeVsDistance_Normed");
            dir_depl->cd();
        }

        TF1 *depl_fit1 = new TF1("depl_fit1","pol1");
        TF1 *depl_fit2 = new TF1("depl_fit2","pol0");

        Float_t *sq_volt = new Float_t[numVolt];
        Float_t *derivative = new Float_t[numVolt];
        for(int i=0;i<numVolt;i++) sq_volt[i] = sqrt(voltages[i]);

        derivative[0]=0;
        derivative[1]=0;
        Float_t max_der=0;
        for(int i=2;i<numVolt;i++){
            derivative[i] = total_charge[i]-total_charge[i-2];
            if(derivative[i]>max_der) max_der=derivative[i];
        }

        int i_start,i_finish,i_plato;
        bool f_start,f_finish;
        f_start = false;
        f_finish = false;
        for(int i=0;i<numVolt;i++) {
            if(!f_start && derivative[i]>0.4*max_der) {i_start=i; f_start=true;}
            if(!f_finish && f_start && derivative[i]<=0.4*max_der) {i_finish=i-1; f_finish=true;}
            if(f_finish && derivative[i]<=0.4*max_der) {i_plato=i; break;}
        }
        Float_t dv = voltages[1]-voltages[0];
        depl_fit1->SetRange(sqrt(i_start*dv),sqrt(i_finish*dv));
        depl_fit2->SetRange(sqrt(i_plato*dv),sqrt(voltages[numVolt]));

        // Plotting the charge
        TGraph *TotalCg=GraphBuilder(numVolt,sq_volt,total_charge,"Sqrt(Voltage) [sqrt(V)]","total charge [.arb]","total charge");
        depl_fit1->SetParameter(0,total_charge[i_start]);
        depl_fit2->SetParameter(0,total_charge[i_plato]);
        TotalCg->Fit("depl_fit1","RQ");
        TotalCg->Fit("depl_fit2","RQ+");

        char depl[100];
        Float_t depl_volt = (depl_fit2->GetParameter(0)-depl_fit1->GetParameter(0))/(depl_fit1->GetParameter(1)-depl_fit2->GetParameter(1));
        sprintf(depl,"U_{depletion} = %.2f V",depl_volt*depl_volt);

        TotalCg->SetTitle(depl);
        TotalCg->Write("DeplVoltage");


        if(config->CH_PhDiode()) {
            dir_depl_normed->cd();

            max_der=0;
            for(int i=2;i<numVolt;i++){
                derivative[i] = total_charge_normed[i]-total_charge_normed[i-2]; // ppor man's derivative; i-1 often doesnt work, try i-2
                if(derivative[i]>max_der) max_der=derivative[i];
            }

            f_start = false;
            f_finish = false;
            for(int i=0;i<numVolt;i++) {
                if(!f_start && derivative[i]>0.4*max_der) {i_start=i; f_start=true;}
                if(!f_finish && f_start && derivative[i]<=0.4*max_der) {i_finish=i-1; f_finish=true;}
                if(f_finish && derivative[i]<=0.4*max_der) {i_plato=i; break;}
            }
            depl_fit1->SetRange(sqrt(i_start*dv),sqrt(i_finish*dv));
            depl_fit2->SetRange(sqrt(i_plato*dv),sqrt(voltages[numVolt]));

            // Plotting the charge
            TGraph *TotalCg1=GraphBuilder(numVolt,sq_volt,total_charge_normed,"Sqrt(Voltage) [sqrt(V)]","total charge [.arb]","total charge");
            TotalCg1->Fit("depl_fit1","RQ");
            TotalCg1->Fit("depl_fit2","RQ+");

            char depl[100];
            depl_volt = (depl_fit2->GetParameter(0)-depl_fit1->GetParameter(0))/(depl_fit1->GetParameter(1)-depl_fit2->GetParameter(1));
            sprintf(depl,"U_{depletion} = %.2f V",depl_volt*depl_volt);

            TotalCg1->SetTitle(depl);
            TotalCg1->Write("DeplVoltage_Normed");

            dir_depl->cd();
        }

        for(int i=0;i<numVolt;i++) {
            delete cc[i];
            delete cc_norm[i];
            delete ph_charge[i];
        }

        delete total_charge;
        delete total_charge_normed;
        delete sq_volt;
        delete derivative;

        return true;
    }

    bool Scanning::DoEdgeVelocity() {

        if(!CheckEdgeVelocity()) {std::cout<<"No data for velocity profile. Skipping..."<<std::endl; return false;}
        TDirectory *dir_vel = f_rootfile->mkdir("Velocity");
        TDirectory *dir_vel_normed;
        if(config->CH_PhDiode()) dir_vel_normed = f_rootfile->mkdir("Velocity_Normed");
        TDirectory *dir_vel_diode;
        if(config->CH_PhDiode()) dir_vel_diode = f_rootfile->mkdir("Velocity_Diode");
        dir_vel->cd();

        Int_t numVolt,numS;
        Float_t Ss,Sc0;

        Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
        Int_t volt_source=(config->VoltSource());            // select
        Int_t scanning_axis=config->ScAxis()-1;         // select scanning axis (0=x,1=y,2=z)
        Float_t *voltages;

        switch(volt_source)
          {
          case 1: numVolt=stct->NU1; voltages=stct->U1.fArray; break; //u1
          case 2: numVolt=stct->NU2; voltages=stct->U2.fArray; break; //u2
          }

        TGraph **cc = new TGraph*[numVolt];
        TGraph **cc_norm; // charge collection graph
        TGraph **field = new TGraph*[numVolt];
        TGraph **velocity_electrons = new TGraph*[numVolt];
        TGraph **velocity_holes = new TGraph*[numVolt];

        TGraph **field_normed = new TGraph*[numVolt];
        TGraph **velocity_electrons_normed = new TGraph*[numVolt];
        TGraph **velocity_holes_normed = new TGraph*[numVolt];

        TGraph **velocity_diode = new TGraph*[numVolt];
        TGraph **field_diode = new TGraph*[numVolt];

        TGraph **ph_charge = new TGraph*[numVolt];

        SwitchAxis(scanning_axis,numS,Ss,Sc0);


        //calculate charge profiles for different voltages
        CalculateCharges(ChNumber,volt_source+2,numVolt,scanning_axis,numS,cc,config->FTlow(),config->FTlow()+config->EV_Time());

        //find integration ranges
        TCTWaveform *wf0;
        switch(volt_source)
          {
          case 1: wf0 = stct->Projection(ChNumber,scanning_axis,0,0,0,numVolt-1,0,numS); break;
          case 2: wf0 = stct->Projection(ChNumber,scanning_axis,0,0,0,0,numVolt-1,numS); break;
          }
        TGraph *charge_max_bias = wf0->CCE(config->FTlow(),config->FThigh());
        delete wf0;
        Double_t left_edge,right_edge;
        FindEdges(charge_max_bias,numS,Ss,left_edge,right_edge);
        std::cout<<std::endl;
        std::cout<<"\tFound detector thickness is: "<<right_edge-left_edge<<" micrometers"<<std::endl;
        std::cout<<"\tIntegrating in the range "<<left_edge<<" - "<<right_edge<<std::endl;
        std::cout<<std::endl;

        Double_t *temp_max_bias_charge = charge_max_bias->GetX();
        Int_t ix1=-1;
        Int_t ix2=-1;
        for(Int_t i=0;i<numS;i++) {
            if(temp_max_bias_charge[i]<=left_edge) ix1=i;
            if(temp_max_bias_charge[i]<=right_edge) ix2=i;
        }
        if(ix1==-1) ix1=0;
        if(ix2==-1) ix2=numS-1;

        //calculate laser charge profiles
        if(config->CH_PhDiode()) CalculateCharges(config->CH_PhDiode()-1,volt_source+2,numVolt,scanning_axis,numS,ph_charge,config->FDLow(),config->FDHigh());

        //calculating the velocity profile asuuming integral(E)dx = Vbias
        Double_t eps = 1e-3;
        Double_t *temp_sensor;
        Double_t *temp_vel_h=new Double_t[numS];
        Double_t *temp_vel_el=new Double_t[numS];
        Double_t *temp_field=new Double_t[numS];
        for(int j=0;j<numVolt;j++) { // FIXME constants seems to prefer smaller values towards low electric fields , why?
            temp_sensor = cc[j]->GetY();

            Double_t a=1;
            if(GraphIntegral(cc[j],left_edge,right_edge)<0) a = -a;
            Double_t dx = 0.5*a;
            Double_t sum = 0;
            while(abs(voltages[j]-sum)>eps) {

                for(int i=0; i<numS; i++) {
                    temp_field[i] = BiSectionMethod(eps,-5e3,1e6,temp_sensor[i],a);
                }
                sum=0;
                for(int i=ix1;i<=ix2;i++) sum+=temp_field[i];
                sum*=1e-4*Ss; // conversion from um to cm
                if(voltages[j]-sum>0) a-=dx;
                else a+=dx;
                dx=dx/2;
            }
            std::cout<<"U = "<<voltages[j]<<" norm const: "<<a<<std::endl;
            for(int i=0;i<numS;i++) {
                temp_field[i] = 1e-4*temp_field[i];
                temp_vel_h[i] = 1e4*temp_field[i]*Mu(temp_field[i],1)/(1+Mu(temp_field[i],1)*1e4*temp_field[i]/config->v_sat());
                temp_vel_el[i] = 1e4*temp_field[i]*Mu(temp_field[i],0)/(1+Mu(temp_field[i],0)*1e4*temp_field[i]/config->v_sat());
            }
            field[j] = GraphBuilder(numS,cc[0]->GetX(),temp_field,"scanning distance [#mum]", "Electric Field, [V/#mum]","Electric Field");
            velocity_holes[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_h,"scanning distance [#mum]", "Velocity [cm/s]","Holes Velocity Profile");
            velocity_electrons[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_el,"scanning distance [#mum]", "Velocity [cm/s]","Electrons Velocity Profile");

        }

        //calculating the velocity profile asuuming integral(E)dx = Vbias with normed data
        Float_t *normcoeff = new Float_t[numVolt];
        if(config->CH_PhDiode()) {
            //calculating the mean charge from the photodetector
            cc_norm=NormedCharge(cc,ph_charge,numVolt);

            Double_t *temp_sensor;

            for(int j=0;j<numVolt;j++) {
                temp_sensor = cc_norm[j]->GetY();
                Double_t a=1;
                if(GraphIntegral(cc_norm[j],left_edge,right_edge)<0) a = -a;
                Double_t dx = 0.5*a;
                Double_t sum = 0;
                while(abs(voltages[j]-sum)>eps) {

                    for(int i=0; i<numS; i++) {
                        temp_field[i] = BiSectionMethod(eps,-5e3,1e6,temp_sensor[i],a);
                    }
                    sum=0;
                    for(int i=ix1;i<=ix2;i++) sum+=temp_field[i];
                    sum*=1e-4*Ss;
                    if(voltages[j]-sum>0) a-=dx;
                    else a+=dx;
                    dx=dx/2;
                }
                if(GraphIntegral(cc_norm[j],left_edge,right_edge)<0) normcoeff[j]=-a;
                else normcoeff[j]=a;
                //std::cout<<"U = "<<voltages[j]<<" norm const: "<<a<<std::endl;
                for(int i=0;i<numS;i++) {
                    temp_field[i] = 1e-4*temp_field[i];
                    temp_vel_h[i] = 1e4*temp_field[i]*Mu(temp_field[i],1)/(1+Mu(temp_field[i],1)*1e4*temp_field[i]/config->v_sat());
                    temp_vel_el[i] = 1e4*temp_field[i]*Mu(temp_field[i],0)/(1+Mu(temp_field[i],0)*1e4*temp_field[i]/config->v_sat());
                }
                field_normed[j] = GraphBuilder(numS,cc[0]->GetX(),temp_field,"scanning distance [#mum]", "Electric Field, [V/#mum]","Electric Field");
                velocity_holes_normed[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_h,"scanning distance [#mum]", "Velocity [cm/s]","Holes Velocity Profile");
                velocity_electrons_normed[j] = GraphBuilder(numS,cc[0]->GetX(),temp_vel_el,"scanning distance [#mum]", "Velocity [cm/s]","Electrons Velocity Profile");

            }
        }

        delete temp_vel_h;
        delete temp_vel_el;
        delete temp_field;


        //calculating the velocity profile using photodiode for estimate of N_e,h, Factor 100 to many :(
        if(config->CH_PhDiode()) {

            Double_t temp_sr_general = 0;
            Double_t *temp_diode;
            Double_t temp_sr = 0;
            for(int j=0;j<numVolt;j++) {
                temp_diode = ph_charge[j]->GetY();
                temp_sr = 0;
                for(int i=0;i<numS;i++) {
                    temp_sr += temp_diode[i];
                }
                temp_sr = temp_sr/(numS);
                temp_sr_general+=temp_sr;
            }
            temp_sr_general=temp_sr_general/(numVolt);

            Double_t ampl = config->ampl(); // FIXME needs accurate update!
            Double_t Res_sensor = config->R_sensor();
            Double_t Res_photo = config->R_diode();
            Double_t Response_photo = config->RespPhoto();
            Double_t Eweight = 1./(config->SampleThickness()*1e-4);
            Double_t E_pair = config->E_pair();
            Double_t diode_multi = config->light_split();

            Double_t *temp_sensor_charge;
            Double_t *temp_diode_charge;
            Double_t *temp_velocity_avg = new Double_t[numS];
            Double_t *temp_field = new Double_t[numS];

            Double_t Neh = 0.624*diode_multi*temp_sr_general/Res_photo/Response_photo/E_pair;
            Neh = Neh*0.02146; // Integral over one strip from energy depostion a.f.o. depth for 80 um strip, and lambda = 11 1/cm
            std::cout<<"Photodiode charge: "<<temp_sr_general<<std::endl;
            std::cout<<"Npairs: "<<1.e7*Neh<<std::endl;

            for(int j=0;j<numVolt;j++) {
                temp_sensor_charge = cc[j]->GetY();
                temp_diode_charge = ph_charge[j]->GetY();
                for(int i=0;i<numS;i++) {
                    temp_velocity_avg[i] = -1e9*0.624*temp_sr_general*temp_sensor_charge[i]/(config->EV_Time()*temp_diode_charge[i]*Eweight*ampl*Res_sensor*Neh);
                    temp_field[i] = 1e-4*temp_velocity_avg[i]/(config->mu0_els()+config->mu0_holes());
                }
                velocity_diode[j] = GraphBuilder(numS,cc[0]->GetX(),temp_velocity_avg,"scanning distance [#mum]", "Velocity [cm/s]","Velocity Profile");
                field_diode[j] = GraphBuilder(numS,cc[0]->GetX(),temp_field,"scanning distance [#mum]", "Electric Field, [V/#mum]","Electric Field Profile");
            }

            delete temp_velocity_avg;
            delete temp_field;

        }

        //write separate charges
        if(config->FSeparateCharges()) {
            dir_vel->cd();
            GraphSeparate(numVolt,velocity_holes,"velocity_profiles_holes","scanning distance, [#mum]","Velocity [cm/s]","Holes Velocity Profile","U = ",voltages);
            GraphSeparate(numVolt,velocity_electrons,"velocity_profiles_electrons","scanning distance, [#mum]","Velocity [cm/s]","Electrons Velocity Profile","U = ",voltages);
            GraphSeparate(numVolt,field,"field_profiles","scanning distance, [#mum]","Electric Field, [V/#mum]","Electric Field Profile","U = ",voltages);
            if(config->CH_PhDiode()) {
                dir_vel_normed->cd();
                GraphSeparate(numVolt,velocity_holes_normed,"velocity_profiles_holes","scanning distance, [#mum]","Velocity [cm/s]","Holes Velocity Profile","U = ",voltages);
                GraphSeparate(numVolt,velocity_electrons_normed,"velocity_profiles_electrons","scanning distance, [#mum]","Velocity [cm/s]","Electrons Velocity Profile","U = ",voltages);
                GraphSeparate(numVolt,field_normed,"field_profiles","scanning distance, [#mum]","Electric Field, [V/#mum]","Electric Field Profile","U = ",voltages);
                dir_vel_diode->cd();
                GraphSeparate(numVolt,ph_charge,"photodiode","scanning distance, [#mum]","Charge,[arb.]","Laser Charge vs Distance","U = ",voltages);
                GraphSeparate(numVolt,velocity_diode,"velocity_profiles","scanning distance, [#mum]","Velocity [cm/s]","Sum Velocity Profile","U = ",voltages);
                GraphSeparate(numVolt,field_diode,"field_profiles","scanning distance, [#mum]","Electric Field, [V/#mum]","Electric Field Profile","U = ",voltages);
                dir_vel->cd();
            }
        }

        //plot graphs for different voltages
        MultiGraphWriter(numVolt,field,"scanning distance [#mum]","Electric Field, [V/#mum]","Electric Field","Electric_Field");
        MultiGraphWriter(numVolt,velocity_holes,"scanning distance [#mum]","Velocity, [cm/s]","Holes Velocity Profile","Vel_Holes");
        MultiGraphWriter(numVolt,velocity_electrons,"scanning distance [#mum]","Velocity, [cm/s]","Electrons Velocity Profile","Vel_Electrons");

        //plot graphs for different voltages with normed data
        if(config->CH_PhDiode()) {
            dir_vel_normed->cd();
            MultiGraphWriter(numVolt,field_normed,"scanning distance [#mum]","Electric Field, [V/#mum]","Electric Field","Electric_Field");
            MultiGraphWriter(numVolt,velocity_holes_normed,"scanning distance [#mum]","Velocity, [cm/s]","Holes Velocity Profile","Vel_Holes");
            MultiGraphWriter(numVolt,velocity_electrons_normed,"scanning distance [#mum]","Velocity, [cm/s]","Electrons Velocity Profile","Vel_Electrons");

            //plots normalization coefficient
            TGraph *coeff = GraphBuilder(numVolt,voltages,normcoeff,"Voltage, [V]","A","Norm Coefficient for Different voltages");
            TF1 *fff = new TF1("log0","[0]*log(x)",10,120);
            coeff->Fit("log0","R");
            coeff->Write("CoeffNorm");

            dir_vel_diode->cd();
            MultiGraphWriter(numVolt,velocity_diode,"scanning distance [#mum]","Velocity [cm/s]","Velocity Profiles","VelocityVsDist_diode");
            MultiGraphWriter(numVolt,field_diode,"scanning distance [#mum]","Electric Field, [V/#mum]","Electric Field","FieldVsDist_diode");
            dir_vel->cd();
        }

        return true;

    }

    // this is independent of the sensor (and hence position), as it only uses the photo diode
    bool Scanning::LaserPowerDrop() {
        f_rootfile->cd();

        Double_t dt;
        if(config->Movements_dt()>0) dt = config->Movements_dt()/60.;
        else dt = 0.001;
        Int_t photo_channel = config->CH_PhDiode();

        Int_t N1=stct->Nx;
        Int_t N2=stct->Ny;
        Int_t N3=stct->Nz;
        Int_t N4=stct->NU1;
        Int_t N5=stct->NU2;
        Int_t numS=N1*N2*N3*N4*N5;

        TH1F **test_graph = new TH1F*[numS];
        Double_t *temp_integral = new Double_t[numS];
        Double_t *xxx = new Double_t[numS];
        Double_t temp_width;
        Int_t tlow;
        Int_t thigh;

        Int_t i=0;
        for(int p=0;p<N5;p++) {
            for(int s=0;s<N4;s++) {
                for(int j=0;j<N3;j++) {
                    for(int k=0;k<N2;k++) {
                        for(int m=0;m<N1;m++) {
                            test_graph[i] = stct->GetHA(photo_channel-1,m,k,j,s,p);
                            temp_width = test_graph[i]->GetBinWidth(1);
                            bool found = false;
                            for(int j=0;j<test_graph[i]->GetNbinsX();j++) {
                                if(!found && test_graph[i]->GetBinContent(j)>5) {tlow = j;found=true;}
                                if(found) { if(test_graph[i]->GetBinContent(j)<0) { thigh = j; break; } }
                            }
                            xxx[i] = i*dt;
                            temp_integral[i] = test_graph[i]->Integral((Int_t)tlow,(Int_t)thigh)*temp_width;
                            i++;
                        }
                    }
                }
            }
        }

        GraphBuilder(numS,xxx,temp_integral,"Time [minutes]","Charge [arb.]","Laser Charge vs Time","ChargeVsTime");

        delete xxx;
        delete temp_integral;
    }

    // this is independent of the sensor (and hence position), as it only uses the photo diode
    bool Scanning::BeamSigma() {

        Double_t dt = config->Movements_dt();
        f_rootfile->cd();
        Int_t photo_channel = config->CH_PhDiode();
        Int_t N1=stct->Nx;
        Int_t N2=stct->Ny;
        Int_t N3=stct->Nz;
        Int_t N4=stct->NU1;
        Int_t N5=stct->NU2;
        Int_t numS=N1*N2*N3*N4*N5;

        TH1F **test_graph = new TH1F*[numS];
        Double_t *temp_integral = new Double_t[numS];
        Double_t *xxx = new Double_t[numS];
        Double_t temp_width;
        Int_t tlow;
        Int_t thigh;
        Double_t ch_min,ch_max;
        Double_t sr = 0;

        Int_t i=0;
        for(int p=0;p<N5;p++) {
            for(int s=0;s<N4;s++) {
                for(int j=0;j<N3;j++) {
                    for(int k=0;k<N2;k++) {
                        for(int m=0;m<N1;m++) {
                            test_graph[i] = stct->GetHA(photo_channel-1,m,k,j,s,p);
                            temp_width = test_graph[i]->GetBinWidth(1);
                            bool found = false;
                            for(int j=0;j<test_graph[i]->GetNbinsX();j++) {
                                if(!found && test_graph[i]->GetBinContent(j)>5) {tlow = j;found=true;}
                                if(found) { if(test_graph[i]->GetBinContent(j)<0) { thigh = j; break; } }
                            }
                            xxx[i] = i;
                            temp_integral[i] = test_graph[i]->Integral((Int_t)tlow,(Int_t)thigh)*temp_width;

                            if(i==0) {
                                ch_min = temp_integral[i];
                                ch_max = temp_integral[i];
                            }
                            if(temp_integral[i]<ch_min) ch_min = temp_integral[i];
                            if(temp_integral[i]>ch_max) ch_max = temp_integral[i];
                            sr+=temp_integral[i];
                            i++;
                        }
                    }
                }
            }
        }
        sr = sr/numS;
        Int_t count_out=0;
        Double_t compare = 2*sr-ch_max;
        for(int i=0;i<numS; i++) {
            if(temp_integral[i]<compare) count_out++;
        }

        TH1F* charge_spread;
        if(count_out>0.1*numS) {
            charge_spread = new TH1F("charge","Charge Distribution",50,ch_min-0.05*ch_min,ch_max+0.05*ch_min);
            for(int i=0;i<numS;i++) charge_spread->Fill(temp_integral[i]);
        }
        else {
            ch_min = compare;
            charge_spread = new TH1F("charge","Charge Distribution",50,ch_min-0.05*ch_min,ch_max+0.05*ch_min);
            for(int i=0;i<numS;i++) {
                if(temp_integral[i]>=compare) charge_spread->Fill(temp_integral[i]);
            }
        }

        TString title;
        title.Form("Laser Charge Distribution in %.1f seconds",dt*numS);
        charge_spread->SetTitle(title.Data());
        charge_spread->GetXaxis()->SetTitle("Charge, [arb.]");
        charge_spread->Draw();
        charge_spread->Write("ChargeDistr");

        delete xxx;
        delete temp_integral;

    }

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

    bool Scanning::CheckFocus() {

        std::cout<<"\t- DO Find Focus: "<< config->DO_focus() <<std::endl;

        Int_t NOpt;
        Int_t NSc;
        std::cout<<"\t\t-- Optical Axis: "<< config->OptAxis() <<std::endl;
        switch(config->OptAxis()) {
        case 1: NOpt = stct->Nx; break;
        case 2: NOpt = stct->Ny; break;
        case 3: NOpt = stct->Nz; break;
        }
        if(NOpt>=1) std::cout<<"\t\t\tOptical axis scan contains "<<NOpt<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tOptical axis contains only "<<NOpt<<" points. Not enough for focusing."<<std::endl;
            return false;
        }
        std::cout<<"\t\t-- Scanning Axis: "<< config->ScAxis() <<std::endl;
        switch(config->ScAxis()) {
        case 1: NSc = stct->Nx; break;
        case 2: NSc = stct->Ny; break;
        case 3: NSc = stct->Nz; break;
        }
        if(NSc>10) std::cout<<"\t\t\tScanning axis contains "<<NSc<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tScanning axis contains only "<<NSc<<" point. Not enough for focusing."<<std::endl;
            return false;
        }

        std::cout<<"FocusSearch Data Test Passed. Processing..."<<std::endl;
        return true;
    }

    bool Scanning::CheckEdgeDepletion() {

        std::cout<<"\t- DO Edge Depletion: "<< config->DO_EdgeDepletion() <<std::endl;

        Int_t NSource;
        Int_t NSc;
        std::cout<<"\t\t-- Voltage Source: "<< config->VoltSource() <<std::endl;
        switch(config->VoltSource()) {
        case 1: NSource = stct->NU1; break;
        case 2: NSource = stct->NU2; break;
        }
        if(NSource>=7) std::cout<<"\t\t\tVoltage scan contains "<<NSource<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tVoltage scan contains only "<<NSource<<" points. Not enough for depletion voltage search."<<std::endl;
            return false;
        }
        std::cout<<"\t\t-- Scanning Axis: "<< config->ScAxis() <<std::endl;
        switch(config->ScAxis()) {
        case 1: NSc = stct->Nx; break;
        case 2: NSc = stct->Ny; break;
        case 3: NSc = stct->Nz; break;
        }
        if(NSc>10) std::cout<<"\t\t\tScanning axis contains "<<NSc<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tScanning axis contains only "<<NSc<<" point. Not enough for depletion voltage search."<<std::endl;
            return false;
        }

        std::cout<<"EdgeDepletionVoltage Data Test Passed. Processing..."<<std::endl;
        return true;
    }

    bool Scanning::CheckEdgeVelocity() {

        std::cout<<"\t- DO Edge Velocity profile: "<< config->DO_EdgeDepletion() <<std::endl;

        Int_t NSource;
        Int_t NSc;
        std::cout<<"\t\t-- Voltage Source: "<< config->VoltSource() <<std::endl;
        switch(config->VoltSource()) {
        case 1: NSource = stct->NU1; break;
        case 2: NSource = stct->NU2; break;
        }
        if(NSource>=1) std::cout<<"\t\t\tVoltage scan contains "<<NSource<<" points. OK"<<std::endl;
        std::cout<<"\t\t-- Scanning Axis: "<< config->ScAxis() <<std::endl;
        switch(config->ScAxis()) {
        case 1: NSc = stct->Nx; break;
        case 2: NSc = stct->Ny; break;
        case 3: NSc = stct->Nz; break;
        }
        if(NSc>10) std::cout<<"\t\t\tScanning axis contains "<<NSc<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tScanning axis contains only "<<NSc<<" point. Not enough for velocity profile."<<std::endl;
            return false;
        }

        std::cout<<"EdgeVelocityProfile Data Test Passed. Processing..."<<std::endl;
        return true;
    }

    // to set number of scanning points, step and x0 for each axis.
    void Scanning::SwitchAxis(Int_t sw, Int_t& nPoints, Float_t& step, Float_t& p0) {
        switch(sw)
          {
          case 0: nPoints=stct->Nx; step=stct->dx; p0=stct->x0;  break; //x-axis
          case 1: nPoints=stct->Ny; step=stct->dy; p0=stct->y0; break; //y-axis
          case 2: nPoints=stct->Nz; step=stct->dz; p0=stct->z0; break; //z-axis
          }
    }

    void Scanning::CalculateCharges(Int_t Channel, Int_t Ax, Int_t numAx, Int_t scanning, Int_t numS, TGraph **charges, Float_t tstart, Float_t tfinish) {
        TCTWaveform **wf = new TCTWaveform*[numAx];
        for(int j=0;j<numAx;j++)
        {
            switch(Ax)
            {
            case 0: wf[j]=stct->Projection(Channel,scanning,j,0,0,0,0,numS); break;
            case 1: wf[j]=stct->Projection(Channel,scanning,0,j,0,0,0,numS); break;
            case 2: wf[j]=stct->Projection(Channel,scanning,0,0,j,0,0,numS); break;
            case 3: wf[j]=stct->Projection(Channel,scanning,0,0,0,j,0,numS); break;
            case 4: wf[j]=stct->Projection(Channel,scanning,0,0,0,0,j,numS); break;
            }

            wf[j]->SetTemperatures(stct->T);
            charges[j]=wf[j]->CCE(tstart,tfinish);   //integrate the charge in time window
            //wf[j]->Clear();
            delete wf[j];
        }
    }

    TGraph** Scanning::NormedCharge(TGraph** sensor, TGraph** photodiode, Int_t numP) {

        TGraph** normed_charge = new TGraph*[numP];
        Int_t numS = sensor[0]->GetN();

        Double_t *temp_sensor;
        Double_t *temp_diode;
        Double_t *temp_normed = new Double_t[numS];
        Double_t temp_sr = 0;
        Double_t temp_sr_general = 0;
        for(int j=0;j<numP;j++) {
            temp_diode = photodiode[j]->GetY();
            temp_sr = 0;
            for(int i=0;i<numS;i++) {
                temp_sr += temp_diode[i];
            }
            temp_sr = temp_sr/numS;
            temp_sr_general+=temp_sr;
        }

        temp_sr_general=temp_sr_general/numP;

        for(int j=0;j<numP;j++) {
            temp_sensor = sensor[j]->GetY();
            temp_diode = photodiode[j]->GetY();
            for(int i=0;i<numS;i++) {
                temp_normed[i] = temp_sr_general*temp_sensor[i]/temp_diode[i];
            }
            normed_charge[j] = new TGraph(numS,sensor[0]->GetX(),temp_normed);
        }
        return normed_charge;

    }

    void Scanning::ChargeCorrelationHist(TGraph** sensor, TGraph** photodetector, Int_t numO) {

        Double_t *temp_phdiode;
        Double_t *temp_detector;
        TH2F *corr_charges = new TH2F("CorrDiode-Sensor","Correlation Diode-Sensor",50,0,1,50,0,1);
        corr_charges->GetXaxis()->SetTitle("PhotoDiode charge [arb.]");
        corr_charges->GetYaxis()->SetTitle("Sensor charge [arb.]");
        corr_charges->SetTitle("Correlation Diode-Sensor");
        Double_t xlow,xup,ylow,yup;
        ylow = sensor[0]->GetHistogram()->GetMinimum();
        yup = sensor[0]->GetHistogram()->GetMaximum();
        xlow = photodetector[0]->GetHistogram()->GetMinimum();
        xup = photodetector[0]->GetHistogram()->GetMaximum();
        for(int i=0;i<numO;i++) {
            temp_detector = sensor[i]->GetY();
            temp_phdiode = photodetector[i]->GetY();
            for(int j=0;j<sensor[0]->GetN();j++) {
                if(ylow>temp_detector[j]) ylow = temp_detector[j];
                if(yup<temp_detector[j]) yup = temp_detector[j];
                if(xlow>temp_phdiode[j]) xlow = temp_phdiode[j];
                if(xup<temp_phdiode[j]) xup = temp_phdiode[j];
            }
        }
        corr_charges->GetXaxis()->SetLimits(0.8*xlow,1.2*xup);
        corr_charges->GetYaxis()->SetLimits(0.8*ylow,1.2*yup);
        for(int i=0;i<numO;i++) {
            temp_detector = sensor[i]->GetY();
            temp_phdiode = photodetector[i]->GetY();
            for(int j=0;j<sensor[0]->GetN();j++) {
                corr_charges->Fill(temp_phdiode[j],temp_detector[j]);

            }
        }
        corr_charges->Write("CorrDiode-Sensor");

    }


    // find position of two edge for a fixed voltage (at least fully depleted) and for one optical distance
    void Scanning::FindEdges(TGraph* gr, Int_t numS, Float_t dx, Double_t& left_edge, Double_t& right_edge) {

        TF1 *ff_left=new TF1("ff_left","[2]*TMath::Erfc((x-[0])/[1])-[3]",0,dx*numS);
        TF1 *ff_right=new TF1("ff_right","[2]*TMath::Erf((x-[0])/[1])-[3]",0,dx*numS);

        double *yy_temp;
        Float_t FWHM = config->FFWHM();
        Float_t *yy = new Float_t[numS];
        Int_t i_max, i_min;
        Float_t max, min, der_max, der_min;
        Float_t FitHeight;

        yy_temp = gr->GetY();
        max = yy_temp[0];
        min = yy_temp[0];
        der_max = yy_temp[1]-yy_temp[0];
        der_min = yy_temp[1]-yy_temp[0];
        yy[0] = 0;
        for(int i=1;i<numS;i++) {
            yy[i-1] = yy_temp[i]-yy_temp[i-1];
            if(yy_temp[i]>max) max=yy_temp[i];
            if(yy_temp[i]<min) min=yy_temp[i];
            if(yy[i-1]>der_max) {
                der_max = yy[i-1];
                i_max = i;
            }
            if(yy[i-1]<der_min) {
                der_min = yy[i-1];
                i_min = i;
            }

        }

        FitHeight = (max-min)/2;
        SetFitParameters(ff_left,i_min*dx,FWHM,FitHeight,FitHeight);
        ff_left->SetRange(0,(i_max+i_min)/2*dx);

        SetFitParameters(ff_right,i_max*dx,FWHM,FitHeight,FitHeight);
        ff_right->SetRange((i_max+i_min)/2*dx,dx*numS);

        gr->Fit("ff_left","NRq");
        gr->Fit("ff_right","NRq+");

        left_edge=ff_left->GetParameter(0);
        right_edge=ff_right->GetParameter(0);

        delete yy;

    }

    // find position of two edge for a fixed voltage (at least fully depleted) and for all optical distances
    void Scanning::FindEdges(TGraph** gr, Int_t numP, Int_t numS, Float_t dx, Float_t* left_pos, Float_t* left_width, Float_t* right_pos, Float_t* right_width) {

        TF1 *ff_left=new TF1("ff_left","[2]*TMath::Erfc((x-[0])/[1])-[3]",0,dx*numS);
        ff_left->SetParName(0,"Left edge");
        ff_left->SetParName(1,"#sigma_left");

        TF1 *ff_right=new TF1("ff_right","[2]*TMath::Erf((x-[0])/[1])-[3]",0,dx*numS);
        ff_right->SetParName(0,"Right edge");
        ff_right->SetParName(1,"#sigma_right");

        double *yy_temp;
        Float_t FWHM = config->FFWHM();
        Float_t *yy = new Float_t[numS];
        Int_t i_max, i_min;
        Float_t max, min, der_max, der_min;
        Float_t FitHeight;

        for(Int_t j=0;j<numP;j++) {

            yy_temp = gr[j]->GetY();
            max = yy_temp[0];
            min = yy_temp[0];
            der_max = yy_temp[1]-yy_temp[0];
            der_min = yy_temp[1]-yy_temp[0];
            yy[0] = 0;
            for(int i=1;i<numS;i++) {
                yy[i-1] = yy_temp[i]-yy_temp[i-1];
                if(yy_temp[i]>max) max=yy_temp[i];
                if(yy_temp[i]<min) min=yy_temp[i];
                if(yy[i-1]>der_max) {
                    der_max = yy[i-1];
                    i_max = i;
                }
                if(yy[i-1]<der_min) {
                    der_min = yy[i-1];
                    i_min = i;
                }

            }

            FitHeight = (max-min)/2;
            SetFitParameters(ff_left,i_min*dx,FWHM,FitHeight,FitHeight);
            ff_left->SetRange(0,(i_max+i_min)/2*dx);

            SetFitParameters(ff_right,i_max*dx,FWHM,FitHeight,FitHeight);
            ff_right->SetRange((i_max+i_min)/2*dx,dx*numS);

            gr[j]->Fit("ff_left","Rq");
            gr[j]->Fit("ff_right","Rq+");
            gr[j]->GetFunction("ff_right")->SetLineColor(kBlue);

            left_width[j]=ff_left->GetParameter(1)*2.35/TMath::Sqrt(2);
            left_pos[j]=ff_left->GetParameter(0);
            right_width[j]=ff_right->GetParameter(1)*2.35/TMath::Sqrt(2);
            right_pos[j]=ff_right->GetParameter(0);

        }

        delete yy;

    }

    TGraph* Scanning::GraphBuilder(Int_t N, Float_t* x, Float_t* y,const char* namex,const char* namey, const char* title) {
        TGraph *temp=new TGraph(N,x,y);
        temp->SetMarkerStyle(21);
        temp->Draw("AP");
        temp->GetHistogram()->GetXaxis()->SetTitle(namex);
        temp->GetHistogram()->GetYaxis()->SetTitle(namey);
        temp->SetTitle(title);
        return temp;
    }

    TGraph* Scanning::GraphBuilder(Int_t N, Float_t* x, Float_t* y,const char* namex,const char* namey, const char* title, const char* write_name) {
        TGraph *temp=new TGraph(N,x,y);
        temp->SetMarkerStyle(21);
        temp->Draw("AP");
        temp->GetHistogram()->GetXaxis()->SetTitle(namex);
        temp->GetHistogram()->GetYaxis()->SetTitle(namey);
        temp->SetTitle(title);
        temp->Write(write_name);
        return temp;
    }

    TGraph* Scanning::GraphBuilder(Int_t N, Double_t* x, Double_t* y,const char* namex,const char* namey, const char* title) {
        TGraph *temp=new TGraph(N,x,y);
        temp->SetMarkerStyle(21);
        temp->Draw("AP");
        temp->GetHistogram()->GetXaxis()->SetTitle(namex);
        temp->GetHistogram()->GetYaxis()->SetTitle(namey);
        temp->SetTitle(title);
        return temp;
    }

    TGraph* Scanning::GraphBuilder(Int_t N, Double_t* x, Double_t* y,const char* namex,const char* namey, const char* title, const char* write_name) {
        TGraph *temp=new TGraph(N,x,y);
        temp->SetMarkerStyle(21);
        temp->Draw("AP");
        temp->GetHistogram()->GetXaxis()->SetTitle(namex);
        temp->GetHistogram()->GetYaxis()->SetTitle(namey);
        temp->SetTitle(title);
        temp->Write(write_name);
        return temp;
    }

    void Scanning::GraphSeparate(Int_t N, TGraph **gr, const char *dir_name, const char *namex, const char *namey, const char *title, const char *name_0, Double_t *name_1) {

        TDirectory* main = TDirectory::CurrentDirectory();
        TDirectory* new_dir = main->mkdir(dir_name);
        new_dir->cd();

        for(int j=0;j<N;j++) {
            std::stringstream ss;
            ss<<name_0<<name_1[j];
            std::string name = ss.str();
            gr[j]->SetTitle(title);
            gr[j]->SetMarkerSize(0);
            gr[j]->SetLineWidth(2.5);
            gr[j]->SetLineColor(kBlack);
            gr[j]->GetXaxis()->SetTitle(namex);
            gr[j]->GetYaxis()->SetTitle(namey);
            gr[j]->Write(name.c_str());
        }

        main->cd();

    }

    void Scanning::GraphSeparate(Int_t N, TGraph **gr, const char *dir_name, const char *namex, const char *namey, const char *title, const char *name_0, Float_t *name_1) {

        TDirectory* main = TDirectory::CurrentDirectory();
        TDirectory* new_dir = main->mkdir(dir_name);
        new_dir->cd();

        for(int j=0;j<N;j++) {
            std::stringstream ss;
            ss<<name_0<<name_1[j];
            std::string name = ss.str();
            gr[j]->SetTitle(title);
            gr[j]->SetMarkerSize(0);
            gr[j]->SetLineWidth(2.5);
            gr[j]->SetLineColor(kBlack);
            gr[j]->GetXaxis()->SetTitle(namex);
            gr[j]->GetYaxis()->SetTitle(namey);
            gr[j]->Write(name.c_str());
        }

        main->cd();

    }

    void Scanning::MultiGraphWriter(Int_t N, TGraph **gr, const char *namex, const char *namey, const char *title, const char *write_name) {

        //TCanvas *canva = new TCanvas(write_name,title,640,480);
        //canva->cd();
        TMultiGraph *mg = new TMultiGraph();

        for(int j=0;j<N;j++)
          {
            gStyle->SetOptFit(0);
            gr[j]->SetLineColor(j%8+1);
            gr[j]->SetMarkerSize(0);
            mg->Add(gr[j],"l");
          }

        mg->Draw("A");
        mg->SetName(write_name);
        mg->SetTitle(title);
        mg->GetXaxis()->SetTitle(namex);
        mg->GetYaxis()->SetTitle(namey);
        mg->Write();
        //canva->Write(write_name);
        //canva->Close();

    }

    void Scanning::SetFitParameters(TF1* ff,Double_t p0,Double_t p1,Double_t p2,Double_t p3) {
        ff->SetParameter(0,p0);
        ff->SetParameter(1,p1);
        ff->SetParameter(2,p2);
        ff->SetParameter(3,p3);
    }

    Double_t Scanning::GraphIntegral(TGraph *gr, Double_t x1, Double_t x2) {
        Int_t Nall = gr->GetN();
        Double_t *x = gr->GetX();
        Double_t *y = gr->GetY();
        Double_t dx = x[1]-x[0];
        Int_t ix1=-1;
        Int_t ix2=-1;
        for(Int_t i=0;i<Nall;i++) {
            if(x[i]<=x1) ix1=i;
            if(x[i]<=x2) ix2=i;
        }
        if(ix1==-1) ix1=0;
        if(ix2==-1) ix2=Nall-1;
        Double_t sum=0;
        for(Int_t i=ix1;i<=ix2;i++) sum+=y[i];
        sum = sum*dx;
        //std::cout<<"ix1 = "<<ix1<<" ix2 = "<<ix2<<std::endl;
        return sum;
    }

    Double_t Scanning::Mu(Double_t E, Int_t Type) {

        if(Type==1) return config->mu0_holes()/(1.+config->mu0_holes()*E/config->v_sat());
        if(Type==0) return config->mu0_els()/sqrt(1.+config->mu0_els()*E/config->v_sat()*config->mu0_els()*E/config->v_sat());

    }

    Double_t Scanning::ff(Double_t E, Double_t Uuu, Double_t a) {
        return Uuu-a*(Mu(E,1)+Mu(E,0))*E;
    }

    Double_t Scanning::BiSectionMethod(Double_t eps, Double_t x1, Double_t x2, Double_t Uuu, Double_t a) {
        if(ff(x1,Uuu,a)==0) return x1;
        if(ff(x2,Uuu,a)==0) return x2;
        Double_t dx = x2-x1;
        Double_t xi;
        while(abs(ff(xi,Uuu,a))>eps) {
            dx = dx/2;
            xi = x1+dx;
            if(sgn(ff(x1,Uuu,a))!=sgn(ff(xi,Uuu,a))) continue;
            else x1=xi;
        }
        //std::cout<<"x = "<<xi<<" Fe = "<<ff(xi,Uuu,a)<<std::endl;
        return xi;
    }

    template <typename T> int Scanning::sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
    Double_t Scanning::abs(Double_t x) {
        if(x<0) return -x;
        return x;
    }
}


