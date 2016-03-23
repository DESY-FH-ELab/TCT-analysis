/**
 * \file
 * \brief Implementation of TCT::ModuleTopMobility methods.
 */

#include "modules/ModuleTopMobility.h"
#include "TStyle.h"

namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleTopMobility::Analysis() {

    Int_t numVolt;

    Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
    Int_t volt_source=(config->VoltSource());            // select
    Float_t *voltages;

    switch(volt_source)
    {
    case 1: numVolt=stct->NU1; voltages=stct->U1.fArray; break; //u1
    case 2: numVolt=stct->NU2; voltages=stct->U2.fArray; break; //u2
    }

    TH1F* ha;

    Float_t *speed = new Float_t[numVolt];
    Float_t *field = new Float_t[numVolt];

    for(int i=0;i<numVolt;i++) {
        ha = stct->GetHA(ChNumber,0,0,0,i,0);
        ha->Draw();

        Int_t tlow = -1;
        Int_t thigh = -1;
        Int_t thalf_low;
        Int_t thalf_high;
        Float_t threshold=5;

        bool found = false;
        Float_t max = ha->GetMaximum();
        Int_t max_bin = ha->GetMaximumBin();
        for(int j=0;j<ha->GetNbinsX();j++) {
            if(!found && ha->GetBinContent(j)>threshold) {tlow = j;found=true;}
            if(j<max_bin && ha->GetBinContent(j)<=0.5*max) thalf_low = j;
            if(j>max_bin && ha->GetBinContent(j)>=0.5*max) thalf_high = j;
            if(found) { if(ha->GetBinContent(j)<0) { thigh = j; break; } }
        }
        if(thigh==-1) thigh = ha->GetNbinsX();
        Float_t temp_width = ha->GetBinWidth(1);
        //std::cout<<"thalflow: "<<ha->GetBinCenter(0)+thalf_low*temp_width<<" thalfhigh: "<<ha->GetBinCenter(0)+thalf_high*temp_width<<std::endl;
        speed[i] = config->SampleThickness()*1e5/(temp_width*(thalf_high-thalf_low));
        //voltages[i] = sqrt(wf1[i]->GetVoltage(0));
        field[i] = voltages[i]/config->SampleThickness()*10000;
        //std::cout<<"Voltage: "<<wf1[i]->GetVoltage(0)<<" Charge: "<<ha->Integral(tlow,thigh)<<" Time: "<<(temp_width*(thigh-tlow))<<std::endl;
    }

    TGraph* velocity_plot = GraphBuilder(numVolt,field,speed,"Electric Field, [V/cm]","Speed, [cm/s]","Charge Carrier Speed vs Electric Field");

    TF1 *mobility_fit = new TF1("mobility_fit","[0]*x/(1.0+[0]*x/[1])");
    mobility_fit->FixParameter(1,config->v_sat());
    mobility_fit->SetParName(0,"Mobility");
    mobility_fit->SetParName(1,"V_sat");
    mobility_fit->SetRange(velocity_plot->GetHistogram()->GetBinCenter(1),velocity_plot->GetHistogram()->GetBinCenter(velocity_plot->GetHistogram()->GetNbinsX()-1));


    velocity_plot->Fit("mobility_fit","RQ");
    gStyle->SetOptFit(1);
    velocity_plot->Write("VelocityVsField");

    delete velocity_plot;
    delete mobility_fit;

    delete speed;
    delete field;
    return true;
}

/// Check before analysis (should be reimplemented by developer)
bool ModuleTopMobility::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

    Int_t NSource;

    std::cout<<"\t\t-- Voltage Source: "<< config->VoltSource() <<std::endl;
    switch(config->VoltSource()) {
    case 1: NSource = stct->NU1; break;
    case 2: NSource = stct->NU2; break;
    }
    if(NSource>=7) std::cout<<"\t\t\tVoltage scan contains "<<NSource<<" points. OK"<<std::endl;
    else {
        std::cout<<"\t\t\tVoltage scan contains only "<<NSource<<" points. Not enough for "<<GetTitle()<<std::endl;
        return false;
    }

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}

}


