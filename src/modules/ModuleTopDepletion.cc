/**
 * \file
 * \brief Implementation of TCT::ModuleTopDepletion methods.
 */

#include "modules/ModuleTopDepletion.h"
#include "TStyle.h"

namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleTopDepletion::Analysis() {

    Int_t numVolt;

    Int_t ChNumber=(config->CH_Det())-1;              // select the oscilloscope channel
    Int_t volt_source=(config->VoltSource());            // select
    Float_t *voltages;

    switch(volt_source)
    {
    case 1: numVolt=stct->NU1; voltages=stct->U1.fArray; break; //u1
    case 2: numVolt=stct->NU2; voltages=stct->U2.fArray; break; //u2
    }

    Float_t* total_charge = new Float_t[numVolt];
    Float_t* total_charge_normed = new Float_t[numVolt];

    TGraph **cc = new TGraph*[1];
    TGraph **cc_norm; // charge collection graph
    TGraph **ph_charge = new TGraph*[1];

    //calculate charge profiles for different voltages
    CalculateCharges(ChNumber,0,1,volt_source+2,numVolt,cc,config->FTlow(),config->FThigh());

    //calculate laser charge profiles
    if(config->CH_PhDiode()) CalculateCharges(config->CH_PhDiode()-1,0,1,volt_source+2,numVolt,ph_charge,config->FDLow(),config->FDHigh());

    //calculating the normed charge distribution
    if(config->CH_PhDiode()) cc_norm=NormedCharge(cc,ph_charge,1);

    //write separate charges

    TF1 *depl_fit1 = new TF1("depl_fit1","pol1");
    TF1 *depl_fit2 = new TF1("depl_fit2","pol0");

    Float_t *sq_volt = new Float_t[numVolt];
    Float_t *derivative = new Float_t[numVolt];

    Double_t *charge_double = cc[0]->GetY();
    Double_t *charge_normed_double;
    if(config->CH_PhDiode()) charge_normed_double = cc_norm[0]->GetY();
    for(int i=0;i<numVolt;i++) {
        total_charge[i] = (Float_t)charge_double[i];
        if(config->CH_PhDiode()) total_charge_normed[i] = (Float_t)charge_normed_double[i];
        sq_volt[i] = sqrt(voltages[i]);
    }

    derivative[0]=0;
    derivative[1]=0;
    Float_t max_der=0;
    Float_t DEPL_THRESHOLD = 0.4;
    for(int i=2;i<numVolt;i++){
        derivative[i] = total_charge[i]-total_charge[i-2];
        if(derivative[i]>max_der) max_der=derivative[i];
    }

    int i_start,i_finish,i_plato;
    bool f_start,f_finish;
    f_start = false;
    f_finish = false;
    for(int i=0;i<numVolt;i++) {
        if(!f_start && derivative[i]>DEPL_THRESHOLD*max_der) {i_start=i; f_start=true;}
        if(!f_finish && f_start && derivative[i]<=DEPL_THRESHOLD*max_der) {i_finish=i-1; f_finish=true;}
        if(f_finish && derivative[i]<=DEPL_THRESHOLD*max_der) {i_plato=i; break;}
    }
    Float_t dv = voltages[1]-voltages[0];
    depl_fit1->SetRange(sqrt(i_start*dv),sqrt(i_finish*dv));
    depl_fit2->SetRange(sqrt(i_plato*dv),sqrt(voltages[numVolt-1]));

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
    delete TotalCg;


    if(config->CH_PhDiode()) {


        max_der=0;
        for(int i=2;i<numVolt;i++){
            derivative[i] = total_charge_normed[i]-total_charge_normed[i-2]; // ppor man's derivative; i-1 often doesnt work, try i-2
            if(derivative[i]>max_der) max_der=derivative[i];
        }

        f_start = false;
        f_finish = false;
        for(int i=0;i<numVolt;i++) {
            if(!f_start && derivative[i]>DEPL_THRESHOLD*max_der) {i_start=i; f_start=true;}
            if(!f_finish && f_start && derivative[i]<=DEPL_THRESHOLD*max_der) {i_finish=i-1; f_finish=true;}
            if(f_finish && derivative[i]<=DEPL_THRESHOLD*max_der) {i_plato=i; break;}
        }
        depl_fit1->SetRange(sqrt(i_start*dv),sqrt(i_finish*dv));
        depl_fit2->SetRange(sqrt(i_plato*dv),sqrt(voltages[numVolt-1]));

        // Plotting the charge
        TGraph *TotalCg1=GraphBuilder(numVolt,sq_volt,total_charge_normed,"Sqrt(Voltage) [sqrt(V)]","total charge [.arb]","total charge");
        TotalCg1->Fit("depl_fit1","RQ");
        TotalCg1->Fit("depl_fit2","RQ+");

        char depl[100];
        depl_volt = (depl_fit2->GetParameter(0)-depl_fit1->GetParameter(0))/(depl_fit1->GetParameter(1)-depl_fit2->GetParameter(1));
        sprintf(depl,"U_{depletion} = %.2f V",depl_volt*depl_volt);

        TotalCg1->SetTitle(depl);
        TotalCg1->Write("DeplVoltage_Normed");
        delete TotalCg1;

    }

    delete cc;
    if(config->CH_PhDiode()) delete cc_norm;
    delete ph_charge;

    delete total_charge;
    delete total_charge_normed;
    delete sq_volt;
    delete derivative;

    return true;
}

/// Check before analysis (should be reimplemented by developer)
bool ModuleTopDepletion::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

    Int_t NSource;

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

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}

}


