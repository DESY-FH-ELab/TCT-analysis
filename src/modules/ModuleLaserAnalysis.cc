/**
 * \file
 * \brief Implementation of TCT::ModuleLaserAnalysis methods.
 */

#include "modules/ModuleLaserAnalysis.h"
#include "TStyle.h"

namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleLaserAnalysis::Analysis() {

    LaserPowerDrop();
    BeamSigma();

    return true;
}

/// Laser charge with time
    bool ModuleLaserAnalysis::LaserPowerDrop() {
        /// This method is independent of the sensor (and hence position), as it only uses the photo diode
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
        delete test_graph;

        return true;
    }

/// Laser Charge spread
    bool ModuleLaserAnalysis::BeamSigma() {
        /// this is independent of the sensor (and hence position), as it only uses the photo diode
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
        delete test_graph;

        return true;
    }

}


