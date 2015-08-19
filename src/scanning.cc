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
#include "TSystem.h"
#include "TFitResultPtr.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TROOT.h"
// External includes
#include "TCTScan.h"

namespace TCT {
    bool Scanning::ReadTCT(char* filename, analysis *ana, bool HasSubs) {
        PSTCT* stct = new PSTCT(filename,-3,2);
        stct->CorrectBaseLine(10.);

        if(!CheckData(stct,ana)) {std::cout<<"File "<<filename<<" contains not enough data for selected operations. Skipping."<<std::endl; return false;}

        std::string outfolder	= ana->OutFolder();
        std::string outpath  = outfolder + "/" + ana->OutSample_ID();
        gSystem->MakeDirectory(outpath.c_str());
        if(HasSubs) outpath  = outpath + "/" + ana->OutSubFolder();
        if(HasSubs) gSystem->MakeDirectory(outpath.c_str());

        std::string pathandfilename;
        if(HasSubs) pathandfilename = outpath  + "/" + ana->OutSample_ID() + "_" + ana->OutSubFolder() + "_" + ana->OutSubsubFolder() + "_";
        else pathandfilename = outpath  + "/" + ana->OutSample_ID() + "_";
        char name[100];
        sprintf(name,"%02d.%02d.%02d_%02d.%02d.%02d",stct->Date[0],stct->Date[1],stct->Date[2],stct->Date[3],stct->Date[4],stct->Date[5]);
        pathandfilename = pathandfilename + name + ".root";

        std::cout << "Output file was created: " << pathandfilename << std::endl;
        TFile* f_rootfile = new TFile(pathandfilename.c_str(),"RECREATE","TCTanalyser");
        f_rootfile->cd();

        if(ana->DO_focus()) DoFocus(f_rootfile,stct,ana);

        f_rootfile->Close();
        stct->~PSTCT();
        delete stct;

        return true;
    }

    bool Scanning::DoFocus(TFile* f_rootfile, PSTCT *stct, analysis* ana) {

        if(!CheckFocus(stct,ana)) {std::cout<<"No data for focusing. Skipping..."<<std::endl; return false;}
        TDirectory *dir_fsearch = f_rootfile->mkdir("FocusSearch");
        dir_fsearch->cd();
        Int_t numO,numS;
        Int_t numWF=1000;

        MeasureWF *wf[numWF];               // waveforms
        MeasureWF *wf1[numWF];

        TGraph *cc[numWF];                  // charge collection graph
        TGraph *cc_norm[numWF];

        TH1F *sample_hist;
        Float_t width[numWF],pos[numWF], abs_pos[numWF];      // array of FWHM values and transition positions
        Float_t width_new[numWF],pos_new[numWF];
        Float_t optical_axis_co[numWF];     // optical axis coordinate
        Float_t minQ[numWF], minQ_new[numWF];
        Float_t strip_w[numWF], strip_w_new[numWF];
        Float_t Ss,Os;                       // Optical axis step
        Float_t Sc0,Opt0;

        TGraph* ph_charge[numWF];

        Int_t ChNumber=(ana->CH_Det())-1;              // select the oscilloscope channel
        Int_t optic_axis=(ana->OptAxis())-1;            // select optic axis (0=x,1=y,2=z)
        Int_t scanning_axis=ana->FPerp()-1;         // select scanning axis (0=x,1=y,2=z)

        f_rootfile->mkdir("sample_signals");
        f_rootfile->cd("sample_signals");

        for(int i=0;i<4;i++) {
            if(stct->WFOnOff[i]) {
                sample_hist = stct->GetHA(i,0,0,0,0,0);
                sample_hist->Write();
            }
        }

        dir_fsearch->cd();

        switch(optic_axis)
          {
          case 0: numO=stct->Nx-1; Os=stct->dx; Opt0=stct->x0; break; //x-axis
          case 1: numO=stct->Ny-1; Os=stct->dy; Opt0=stct->y0; break; //y-axis
          case 2: numO=stct->Nz-1; Os=stct->dz; Opt0=stct->z0; break; //z-axis
          }
        switch(scanning_axis)
          {
          case 0: numS=stct->Nx-1; Ss=stct->dx; Sc0=stct->x0;  break; //x-axis
          case 1: numS=stct->Ny-1; Ss=stct->dx; Sc0=stct->y0; break; //y-axis
          case 2: numS=stct->Nz-1; Ss=stct->dx; Sc0=stct->z0; break; //z-axis
          }
        for(int j=0;j<=numO;j++)
           {
             switch(optic_axis)
           {
           case 0: wf[j]=stct->Projection(ChNumber,scanning_axis,j,0,0,0,0,numS+1); break;
           case 1: wf[j]=stct->Projection(ChNumber,scanning_axis,0,j,0,0,0,numS+1); break;
           case 2: wf[j]=stct->Projection(ChNumber,scanning_axis,0,0,j,0,0,numS+1); break;
           }

             wf[j]->SetTemperatures(stct->T);
             wf[j]->ScaleHisto(0.0001); //use any value to set the scale
             //wf[j]->CorrectBaseLine(1); //correct base line - not needed
             cc[j]=wf[j]->CCE(ana->FTlow(),ana->FThigh());   //integrate the charge in time window 0-80 ns

           }

        //laser charge distribution
        if(ana->CH_PhDiode()) {


            for(int j=0;j<=numO;j++)
            {
                switch(optic_axis)
                {
                case 0: wf1[j]=stct->Projection(ana->CH_PhDiode()-1,scanning_axis,j,0,0,0,0,numS+1); break;
                case 1: wf1[j]=stct->Projection(ana->CH_PhDiode()-1,scanning_axis,0,j,0,0,0,numS+1); break;
                case 2: wf1[j]=stct->Projection(ana->CH_PhDiode()-1,scanning_axis,0,0,j,0,0,numS+1); break;
                }

                wf1[j]->SetTemperatures(stct->T);
                //wf1[j]->ScaleHisto(0.0001); //use any value to set the scale
                //wf[j]->CorrectBaseLine(1); //correct base line - not needed
                ph_charge[j]=wf1[j]->CCE(ana->FDLow(),ana->FDHigh());   //integrate the charge in time window 0-80 ns
            }

            Double_t *temp_phdiode;
            Double_t *temp_detector;
            TH2F *corr_charges = new TH2F("CorrDiode-Sensor","Correlation Diode-Sensor",50,0,1,50,0,1);
            corr_charges->GetXaxis()->SetTitle("PhotoDiode charge [arb.]");
            corr_charges->GetYaxis()->SetTitle("Sensor charge [arb.]");
            corr_charges->SetTitle("Correlation Diode-Sensor");
            Double_t xlow,xup,ylow,yup;
            ylow = cc[0]->GetHistogram()->GetMinimum();
            yup = cc[0]->GetHistogram()->GetMaximum();
            xlow = ph_charge[0]->GetHistogram()->GetMinimum();
            xup = ph_charge[0]->GetHistogram()->GetMaximum();
            for(int i=0;i<=numO;i++) {
                temp_detector = cc[i]->GetY();
                temp_phdiode = ph_charge[i]->GetY();
                for(int j=0;j<=numS;j++) {
                    if(ylow>temp_detector[j]) ylow = temp_detector[j];
                    if(yup<temp_detector[j]) yup = temp_detector[j];
                    if(xlow>temp_phdiode[j]) xlow = temp_phdiode[j];
                    if(xup<temp_phdiode[j]) xup = temp_phdiode[j];
                }
            }
            corr_charges->GetXaxis()->SetLimits(0.8*xlow,1.2*xup);
            corr_charges->GetYaxis()->SetLimits(0.8*ylow,1.2*yup);
            for(int i=0;i<=numO;i++) {
                temp_detector = cc[i]->GetY();
                temp_phdiode = ph_charge[i]->GetY();
                for(int j=0;j<=numS;j++) {
                    corr_charges->Fill(temp_phdiode[j],temp_detector[j]);
                    //if(i<5 && j<3) std::cout<<temp_phdiode[j]<<" "<<temp_detector[j]<<std::endl;
                }
                //if(i<5)std::cout<<std::endl;
            }
            corr_charges->Write("CorrDiode-Sensor");
        }

        // fit definition

        Float_t LowLim=ana->FFitLow();             //   low limit for fit
        Float_t HiLim;
        if((ana->FFitHigh())==-1) HiLim=Sc0+Ss*numS;
        else HiLim=ana->FFitHigh();             //   high limit for fit
        Float_t FWHM=ana->FFWHM();               //   expected FWHM
        Float_t Level=ana->FLevel();               //  start of flat level



        TF1 *ff2=new TF1("ff2","[2]/2.*(TMath::Erfc((x-[0])/[1]) + (TMath::Erf((x-[0]-[3])/[1]) + 1))",LowLim,HiLim);
        ff2->SetParameter(0,Level);
        ff2->SetParameter(1,FWHM);
        ff2->SetParameter(2,cc[0]->GetHistogram()->GetMaximum());
        ff2->SetParameter(3,0.7*FWHM);

        for(int j=0;j<=numO;j++){
            cc[j]->Fit("ff2","Rq");
            gStyle->SetOptFit(1);

            width[j]=ff2->GetParameter(1)*2.35/TMath::Sqrt(2);
            pos[j]=ff2->GetParameter(0)+ ff2->GetParameter(3)/2;
            minQ[j] = ff2->Eval(ff2->GetParameter(0) + ff2->GetParameter(3)/2)/ff2->GetParameter(2);
            strip_w[j] = ff2->GetParameter(3);
            optical_axis_co[j]=Opt0+j*Os;
        }


        //calculating the normed charge distribution
        if(ana->CH_PhDiode()) {

            ff2->SetParameter(0,Level);
            ff2->SetParameter(1,FWHM);
            ff2->SetParameter(2,3);
            ff2->SetParameter(3,0.7*FWHM);

            Double_t *temp_Y1;
            Double_t *temp_Y2;
            Double_t *temp_Y3 = new Double_t[numS+1];
            Double_t temp_sr = 0;
            Double_t temp_sr_general = 0;
            for(int j=0;j<=numO;j++) {
                temp_Y2 = ph_charge[j]->GetY();
                temp_sr = 0;
                for(int i=0;i<=numS;i++) {
                    temp_sr += temp_Y2[i];
                }
                temp_sr = temp_sr/(numS+1);
                temp_sr_general+=temp_sr;
                //std::cout<<" "<<temp_sr;
            }
            temp_sr_general=temp_sr_general/(numO+1);
            //std::cout<<"temp sr general: "<<temp_sr_general<<std::endl;
            for(int j=0;j<=numO;j++) {
                temp_Y1 = cc[j]->GetY();
                temp_Y2 = ph_charge[j]->GetY();
                for(int i=0;i<=numS;i++) {
                    temp_Y3[i] = temp_sr_general*temp_Y1[i]/temp_Y2[i];
                }
                cc_norm[j] = new TGraph(numS+1,cc[0]->GetX(),temp_Y3);
                if(j==0) ff2->SetParameter(2,cc_norm[j]->GetHistogram()->GetMaximum());
                cc_norm[j]->Fit("ff2","Rq");
                gStyle->SetOptFit(1);

                width_new[j]=ff2->GetParameter(1)*2.35/TMath::Sqrt(2);
                pos_new[j]=ff2->GetParameter(0)+ ff2->GetParameter(3)/2;
                minQ_new[j] = ff2->Eval(ff2->GetParameter(0) + ff2->GetParameter(3)/2)/ff2->GetParameter(2);
                strip_w_new[j] = ff2->GetParameter(3);
            }
        }


        if(ana->FSeparateCharges()) {
            TDirectory* dir_charges = dir_fsearch->mkdir("charges");
            dir_charges->cd();

            for(int j=0;j<=numO;j++) {
                std::stringstream ss;
                ss<<"Z = "<<stct->z0+j*stct->dz;
                std::string name = ss.str();
                cc[j]->SetTitle("ChargeVsDist");
                cc[j]->GetXaxis()->SetTitle("scanning distance [#mum]");
                cc[j]->Write(name.c_str());;
            }

            dir_fsearch->cd();

            if(ana->CH_PhDiode()) {
                TDirectory* dir_phdiode = dir_fsearch->mkdir("photodiode");
                dir_phdiode->cd();

                for(int j=0;j<=numO;j++) {
                    std::stringstream ss;
                    ss<<"Z = "<<stct->z0+j*stct->dz;
                    std::string name = ss.str();
                    ph_charge[j]->SetTitle("LaserChargeVsDist");
                    ph_charge[j]->GetXaxis()->SetTitle("scanning distance [#mum]");
                    ph_charge[j]->Write(name.c_str());
                }

                dir_fsearch->cd();

                TDirectory* dir_charges_normed = dir_fsearch->mkdir("charges_normed");
                dir_charges_normed->cd();

                for(int j=0;j<=numO;j++) {
                    std::stringstream ss;
                    ss<<"Z = "<<stct->z0+j*stct->dz;
                    std::string name = ss.str();
                    cc_norm[j]->SetTitle("ChargeVsDist");
                    cc_norm[j]->GetXaxis()->SetTitle("scanning distance [#mum]");
                    cc_norm[j]->Write(name.c_str());
                }

                dir_fsearch->cd();
            }
        }

        //plot graphs at different positions along optical axis
        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("scanning distance [#mum]");

        for(int j=0;j<=numO;j++)
          {
            gStyle->SetOptFit(0);
            cc[j]->SetLineColor(j%8+1);
            if(j==0)
          {
            cc[j]->Draw("AL");
          }
            else cc[j]->Draw("L");
            mg->Add(cc[j]);
          }

        mg->Draw("AP");
        mg->SetName("ChargeVsDist");
        mg->GetYaxis()->SetTitle("Charge [arb.]");
        mg->GetXaxis()->SetTitle("scanning distance [#mum]");
        mg->Write();

        //plot graphs at different positions along optical axis with normed data
        if(ana->CH_PhDiode()) {

            TMultiGraph *mg1 = new TMultiGraph();
            mg1->SetTitle("scanning distance [#mum]");

            for(int j=0;j<=numO;j++)
              {
                cc_norm[j]->SetLineColor(j%8+1);
                if(j==0)
              {
                cc_norm[j]->Draw("AL");
              }
                else cc_norm[j]->Draw("L");
                mg1->Add(cc_norm[j]);
              }

            mg1->Draw("AP");
            mg1->SetName("ChargeVsDist_Normed");
            mg1->GetYaxis()->SetTitle("Charge [arb]");
            mg1->GetXaxis()->SetTitle("scanning distance [#mum]");
            mg1->Write();

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
        TGraph *FWHMg=new TGraph(numO+1,optical_axis_co,width);
        FWHMg->SetMarkerStyle(21);
        FWHMg->Draw("AP");
        FWHMg->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        FWHMg->GetHistogram()->GetYaxis()->SetTitle("FWHM [#mum]");
        ff_pol1->SetParameter(0,FWHMg->GetMean());
        ff_pol1->SetParameter(2,FWHMg->GetMinimum());
        FWHMg->Fit("my_pol1","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg->SetTitle("Gaussian Beam Profile");
        FWHMg->Write("FWHM");


        // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance

        TGraph *MINQg=new TGraph(numO+1,optical_axis_co,minQ);
        MINQg->SetMarkerStyle(21);
        MINQg->Draw("AP");
        MINQg->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        MINQg->GetHistogram()->GetYaxis()->SetTitle("min Q[rel. to max]");
        ff_pol3->SetParameter(0,MINQg->GetMean());
        ff_pol3->SetParameter(2,MINQg->GetMinimum());
        MINQg->Fit("my_pol3","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)MINQg->FindObject("stats");
        st->SetFitFormat(".5g");
        MINQg->SetTitle("Minimum Charge");
        MINQg->Write("MinCharge");

        // Find the missalignment between z and optical axis
        for(int j=0;j<=numO;j++) abs_pos[j]=pos[j]+Sc0;

        TGraph *POSg=new TGraph(numO+1,optical_axis_co,abs_pos);
        POSg->SetMarkerStyle(21);
        POSg->Draw("AP");
        POSg->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        POSg->GetHistogram()->GetYaxis()->SetTitle("position of the edge [#mum]");
        POSg->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg->SetTitle("Missalignment");
        POSg->Write("Missalignment");

        // Plotting best strip width

        TGraph *StrWg=new TGraph(numO+1,optical_axis_co,strip_w);
        StrWg->SetMarkerStyle(21);
        StrWg->Draw("AP");
        StrWg->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        StrWg->GetHistogram()->GetYaxis()->SetTitle("strip width [#mum]");
        //ff_pol3->SetParameter(0,MINQg->GetMean());
       // ff_pol3->SetParameter(2,MINQg->GetMinimum());
        //MINQg->Fit("my_pol3","q");
        //gStyle->SetOptFit(1);
        //st = (TPaveStats*)MINQg->FindObject("stats");
        //st->SetFitFormat(".5g");
        StrWg->SetTitle("Strip Width");
        StrWg->Write("StripWidth");

        if(ana->CH_PhDiode()) {

            //draw the gaussian beam profile with normed charge
            TGraph *FWHMg1=new TGraph(numO+1,optical_axis_co,width_new);
            FWHMg1->SetMarkerStyle(21);
            FWHMg1->Draw("AP");
            FWHMg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
            FWHMg1->GetHistogram()->GetYaxis()->SetTitle("FWHM [#mum]");
            ff_pol2->SetParameter(0,FWHMg1->GetMean());
            ff_pol2->SetParameter(2,FWHMg1->GetMinimum());
            FWHMg1->Fit("my_pol2","q");
            gStyle->SetOptFit(1);
            st = (TPaveStats*)FWHMg1->FindObject("stats");
            st->SetFitFormat(".5g");
            FWHMg1->SetTitle("Gaussian Beam Profile");
            FWHMg1->Write("FWHM_Normed");

            // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance with normed data

            TGraph *MINQg1=new TGraph(numO+1,optical_axis_co,minQ_new);
            MINQg1->SetMarkerStyle(21);
            MINQg1->Draw("AP");
            MINQg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
            MINQg1->GetHistogram()->GetYaxis()->SetTitle("min Q [rel. to max]");
            ff_pol4->SetParameter(0,MINQg1->GetMean());
            ff_pol4->SetParameter(2,MINQg1->GetMinimum());
            MINQg1->Fit("my_pol4","q");
            gStyle->SetOptFit(1);
            st = (TPaveStats*)MINQg1->FindObject("stats");
            st->SetFitFormat(".5g");
            MINQg1->SetTitle("Minimum Charge");
            MINQg1->Write("MinCharge_Normed");

            // Find the missalignment between z and optical axis
            for(int j=0;j<=numO;j++) abs_pos[j]=pos_new[j]+Sc0;

            TGraph *POSg1=new TGraph(numO+1,optical_axis_co,abs_pos);
            POSg1->SetMarkerStyle(21);
            POSg1->Draw("AP");
            POSg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
            POSg1->GetHistogram()->GetYaxis()->SetTitle("position of the edge [#mum]");
            POSg1->Fit("pol1","q");
            gStyle->SetOptFit(1);
            POSg1->SetTitle("Missalignment");
            POSg1->Write("Missalignment_Normed");

            // Plotting best strip width

            TGraph *StrWg1=new TGraph(numO+1,optical_axis_co,strip_w_new);
            StrWg1->SetMarkerStyle(21);
            StrWg1->Draw("AP");
            StrWg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
            StrWg1->GetHistogram()->GetYaxis()->SetTitle("strip width [#mum]");
            //ff_pol3->SetParameter(0,MINQg->GetMean());
           // ff_pol3->SetParameter(2,MINQg->GetMinimum());
            //MINQg->Fit("my_pol3","q");
            //gStyle->SetOptFit(1);
            //st = (TPaveStats*)MINQg->FindObject("stats");
            //st->SetFitFormat(".5g");
            StrWg1->SetTitle("Strip Width");
            StrWg1->Write("StripWidth_Normed");

        }
/*
        std::cout<<"(FWHM)Focus is at: "<<ff_pol1->GetParameter(0)<<std::endl;
        if(ana->CH_PhDiode()) std::cout<<"(FWHM Normed)Focus is at: "<<ff_pol2->GetParameter(0)<<std::endl;
        std::cout<<"(QMin)Focus is at: "<<ff_pol3->GetParameter(0)<<std::endl;
        if(ana->CH_PhDiode()) std::cout<<"(Qmin Normed)Focus is at: "<<ff_pol4->GetParameter(0)<<std::endl;
*/


        // Correlation between Qmin and FWHM
        TGraph *CORg=new TGraph(numO+1,width_new,minQ_new);
        CORg->SetMarkerStyle(21);
        CORg->SetLineWidth(0);
        CORg->Draw("A");
        CORg->GetHistogram()->GetXaxis()->SetTitle("FWHM [#mum]");
        CORg->GetHistogram()->GetYaxis()->SetTitle("min Q from fit [arb]");
        //CORg->Fit("pol2");
        CORg->SetTitle("Correlation between FWHM and Qmin");
        CORg->Write("Corr_FWHM_Q");


        return true;
    }

    bool Scanning::SimulateDoFocus(TFile* f_rootfile, analysis* ana) {

        TDirectory *dir_fsearch = f_rootfile->mkdir("FocusSearch");
        dir_fsearch->cd();
        Int_t numO,numS;
        Int_t numWF=1000;

        TGraph *cc[numWF];                  // charge collection graph
        TGraph *cc_norm[numWF];

        Float_t abs_pos[numWF];      // array of FWHM values and transition positions
        Float_t width_new[numWF],pos_new[numWF];
        Float_t optical_axis_co[numWF];     // optical axis coordinate
        Float_t scanning_axis_co[numWF];     // optical axis coordinate
        Float_t minQ_new[numWF];
        Float_t strip_w_new[numWF];
        Float_t Ss,Os;                       // Optical axis step
        Float_t Sc0,Opt0;
        Float_t Pi = 3.1415926535;

        Opt0 = 8600;
        Os = 50;
        numO = 24;
        Sc0 = 0;
        Ss = 2;
        numS = 40;
        Float_t strip_w = 6.8;
        Float_t data_strip_x = 40;

        Float_t beam_sigma[numWF];
        Float_t maximal_beam_sigma = 15;
        Float_t minimal_beam_sigma = 5.7;
        Float_t norm=0;
        for(int i=1;i<=numO/2;i++) norm+=i*i;
        norm = 2*norm;
        for(int j=0;j<=numO;j++) {
            optical_axis_co[j] = Opt0 + Os*j;
            beam_sigma[j] = minimal_beam_sigma + 6.85e-5/numO*10*(j-numO/2)*(j-numO/2)*Os*Os;
            //std::cout<<beam_sigma[j]<<std::endl;
        }
        for(int j=0;j<=numS;j++) {
            scanning_axis_co[j] = Sc0 + Ss*numS;
        }


        // fit definition

        Float_t LowLim=Sc0;             //   low limit for fit
        Float_t HiLim=Sc0+Ss*numS;
        Float_t FWHM=maximal_beam_sigma/2;               //   expected FWHM
        Float_t Level=data_strip_x;               //  start of flat level



        TF1 *ff2=new TF1("ff2","[2]/2.*(TMath::Erfc((x-[0])/[1]) + (TMath::Erf((x-[0]-[3])/[1]) + 1))",LowLim,HiLim);



        ff2->SetParameter(0,Level);
        ff2->SetParameter(1,FWHM);
        ff2->SetParameter(2,3);
        ff2->SetParameter(3,0.7*FWHM);

        Float_t *temp_Y3 = new Float_t[numS+1];
        Float_t beam_x;

        TF1 *fgaus = new TF1("ff1","gaus",LowLim,HiLim);

        for(int j=0;j<=numO;j++) {
            //while(beam_x<=Sc0+Ss*numS) {
            for(Int_t counts=0;counts<=numS;counts++) {
                beam_x = Sc0+Ss*counts;
                fgaus->SetParameters(1.0/sqrt(2*Pi)/beam_sigma[j],beam_x,beam_sigma[j]);
                temp_Y3[counts] = fgaus->Integral(Sc0-4*beam_sigma[j],data_strip_x-strip_w/2)+fgaus->Integral(data_strip_x+strip_w/2,Sc0+Ss*numS+4*beam_sigma[j]);
                scanning_axis_co[counts] = beam_x;

            }
            cc_norm[j] = new TGraph(numS+1,scanning_axis_co,temp_Y3);
            if(j==0) ff2->SetParameter(2,cc_norm[j]->GetHistogram()->GetMaximum());
            cc_norm[j]->Fit("ff2","Rq");
            gStyle->SetOptFit(1);

            width_new[j]=ff2->GetParameter(1)*2.35/TMath::Sqrt(2);
            pos_new[j]=ff2->GetParameter(0)+ ff2->GetParameter(3)/2;
            minQ_new[j] = ff2->Eval(ff2->GetParameter(0) + ff2->GetParameter(3)/2)/ff2->GetParameter(2);
            strip_w_new[j] = ff2->GetParameter(3);
        }


        if(ana->FSeparateCharges()) {

            TDirectory* dir_charges_normed = dir_fsearch->mkdir("charges_normed");
            dir_charges_normed->cd();

            for(int j=0;j<=numO;j++) {
                std::stringstream ss;
                ss<<"Z = "<<optical_axis_co[j];
                std::string name = ss.str();
                cc_norm[j]->SetTitle("ChargeVsDist");
                cc_norm[j]->GetXaxis()->SetTitle("scanning distance [#mum]");
                cc_norm[j]->Write(name.c_str());
            }

            dir_fsearch->cd();
        }

        //plot graphs at different positions along optical axis with normed data

        TMultiGraph *mg1 = new TMultiGraph();
        mg1->SetTitle("scanning distance [#mum]");

        for(int j=0;j<=numO;j++)
        {
            cc_norm[j]->SetLineColor(j%8+1);
            if(j==0)
            {
                cc_norm[j]->Draw("AL");
            }
            else cc_norm[j]->Draw("L");
            mg1->Add(cc_norm[j]);
        }

        mg1->Draw("AP");
        mg1->SetName("ChargeVsDist_Normed");
        mg1->GetYaxis()->SetTitle("Charge [arb]");
        mg1->GetXaxis()->SetTitle("scanning distance [#mum]");
        mg1->Write();


        TF1 *ff_pol2 = new TF1("my_pol2","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol2->SetParName(0,"Focus Position");
        TF1 *ff_pol4 = new TF1("my_pol4","[1]*(x-[0])*(x-[0])+[2]");
        ff_pol4->SetParName(0,"Focus Position");
        TPaveStats *st;


        //draw the gaussian beam profile with normed charge
        TGraph *FWHMg1=new TGraph(numO+1,optical_axis_co,width_new);
        FWHMg1->SetMarkerStyle(21);
        FWHMg1->Draw("AP");
        FWHMg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        FWHMg1->GetHistogram()->GetYaxis()->SetTitle("FWHM [#mum]");
        ff_pol2->SetParameter(0,FWHMg1->GetMean());
        ff_pol2->SetParameter(2,FWHMg1->GetMinimum());
        FWHMg1->Fit("my_pol2","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg1->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg1->SetTitle("Gaussian Beam Profile");
        FWHMg1->Write("FWHM_Normed");

        // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance with normed data

        TGraph *MINQg1=new TGraph(numO+1,optical_axis_co,minQ_new);
        MINQg1->SetMarkerStyle(21);
        MINQg1->Draw("AP");
        MINQg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        MINQg1->GetHistogram()->GetYaxis()->SetTitle("min Q [rel. to max]");
        ff_pol4->SetParameter(0,MINQg1->GetMean());
        ff_pol4->SetParameter(2,MINQg1->GetMinimum());
        MINQg1->Fit("my_pol4","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)MINQg1->FindObject("stats");
        st->SetFitFormat(".5g");
        MINQg1->SetTitle("Minimum Charge");
        MINQg1->Write("MinCharge_Normed");

        // Find the missalignment between z and optical axis
        for(int j=0;j<=numO;j++) abs_pos[j]=pos_new[j]+Sc0;

        TGraph *POSg1=new TGraph(numO+1,optical_axis_co,abs_pos);
        POSg1->SetMarkerStyle(21);
        POSg1->Draw("AP");
        POSg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        POSg1->GetHistogram()->GetYaxis()->SetTitle("position of the edge [#mum]");
        POSg1->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg1->SetTitle("Missalignment");
        POSg1->Write("Missalignment_Normed");

        // Plotting best strip width

        TGraph *StrWg1=new TGraph(numO+1,optical_axis_co,strip_w_new);
        StrWg1->SetMarkerStyle(21);
        StrWg1->Draw("AP");
        StrWg1->GetHistogram()->GetXaxis()->SetTitle("optical distance [#mum]");
        StrWg1->GetHistogram()->GetYaxis()->SetTitle("strip width [#mum]");
        //ff_pol3->SetParameter(0,MINQg->GetMean());
        // ff_pol3->SetParameter(2,MINQg->GetMinimum());
        //MINQg->Fit("my_pol3","q");
        //gStyle->SetOptFit(1);
        //st = (TPaveStats*)MINQg->FindObject("stats");
        //st->SetFitFormat(".5g");
        StrWg1->SetTitle("Strip Width");
        StrWg1->Write("StripWidth_Normed");

/*
        std::cout<<"(FWHM)Focus is at: "<<ff_pol1->GetParameter(0)<<std::endl;
        if(ana->CH_PhDiode()) std::cout<<"(FWHM Normed)Focus is at: "<<ff_pol2->GetParameter(0)<<std::endl;
        std::cout<<"(QMin)Focus is at: "<<ff_pol3->GetParameter(0)<<std::endl;
        if(ana->CH_PhDiode()) std::cout<<"(Qmin Normed)Focus is at: "<<ff_pol4->GetParameter(0)<<std::endl;
*/


        // Correlation between Qmin and FWHM
        TGraph *CORg=new TGraph(numO+1,width_new,minQ_new);
        CORg->SetMarkerStyle(21);
        CORg->Draw("P");
        CORg->GetHistogram()->GetXaxis()->SetTitle("FWHM [#mum]");
        CORg->GetHistogram()->GetYaxis()->SetTitle("min Q from fit [arb]");
        //CORg->Fit("pol2");
        CORg->SetTitle("Correlation between FWHM and Qmin");
        CORg->Write("Corr_FWHM_Q");


        return true;
    }

    bool Scanning::CheckData(PSTCT *stct1, analysis *ana) {

        std::cout<<"Checking Channels set:"<<std::endl;
        std::cout<<"\t- Detector Signal Channel: "<< ana->CH_Det() <<std::endl;
        if(ana->CH_Det()) {
            if(stct1->WFOnOff[ana->CH_Det()-1]) std::cout<<"\t\t Data OK"<<std::endl;
            else {
                std::cout<<"\t\t This channel has no data. Please check the settings."<<std::endl;
                return false;
            }
        }
        std::cout<<"\t- Trigger Channel: "<< ana->CH_Trig() <<std::endl;
        if(ana->CH_Trig()) {
            if(stct1->WFOnOff[ana->CH_Trig()-1]) std::cout<<"\t\t Data OK"<<std::endl;
            else {
                std::cout<<"\t\t This channel has no data. Please check the settings."<<std::endl;
                return false;
            }
        }
        std::cout<<"\t- Photodiode Channel: "<< ana->CH_PhDiode() <<std::endl;
        if(ana->CH_PhDiode()) {
            if(stct1->WFOnOff[ana->CH_PhDiode()-1]) std::cout<<"\t\t Data OK"<<std::endl;
            else {
                std::cout<<"\t\t This channel has no data. Please check the settings."<<std::endl;
                return false;
            }
        }

        std::cout<<"Channel test passed. Processing..."<<std::endl;
        return true;
    }

    bool Scanning::CheckFocus(PSTCT *stct1, analysis* ana) {

        std::cout<<"\t- DO Find Focus: "<< ana->DO_focus() <<std::endl;

        Int_t NOpt;
        Int_t NSc;
        std::cout<<"\t\t-- Optical Axis: "<< ana->OptAxis() <<std::endl;
        switch(ana->OptAxis()) {
        case 1: NOpt = stct1->Nx; break;
        case 2: NOpt = stct1->Ny; break;
        case 3: NOpt = stct1->Nz; break;
        }
        if(NOpt>=5) std::cout<<"\t\t\tOptical axis scan contains "<<NOpt<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tOptical axis contains only "<<NOpt<<" points. Not enough for focusing."<<std::endl;
            return false;
        }
        std::cout<<"\t\t-- Scanning Axis: "<< ana->FPerp() <<std::endl;
        switch(ana->FPerp()) {
        case 1: NSc = stct1->Nx; break;
        case 2: NSc = stct1->Ny; break;
        case 3: NSc = stct1->Nz; break;
        }
        if(NSc>1) std::cout<<"\t\t\tScanning axis contains "<<NSc<<" points. OK"<<std::endl;
        else {
            std::cout<<"\t\t\tScanning axis contains only "<<NSc<<" point. Not enough for focusing."<<std::endl;
            return false;
        }

        std::cout<<"FocusSearch Data Test Passed. Processing..."<<std::endl;
        return true;
    }
}
