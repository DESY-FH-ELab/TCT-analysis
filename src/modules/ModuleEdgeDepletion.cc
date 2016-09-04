/**
 * \file
 * \brief Implementation of TCT::ModuleEdgeDepletion methods.
 */

#include "modules/ModuleEdgeDepletion.h"
#include "TStyle.h"

namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleEdgeDepletion::Analysis() {

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

    Float_t* total_charge = new Float_t[numVolt];
    Float_t* total_charge_normed = new Float_t[numVolt];

    TGraph **cc = new TGraph*[numVolt];
    TGraph **cc_norm; // charge collection graph
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


    //std::cout<<"charge left: "<<charge_left<<" charge right: "<<charge_right<<std::endl;
    for(int j=0;j<numVolt;j++) {
        total_charge[j] = GraphIntegral(cc[j],left_edge,right_edge);
    }
    if(total_charge[numVolt-1]<0) {
        for(int j=0;j<numVolt;j++) total_charge[j]=-total_charge[j];
    }
    if(config->CH_PhDiode()) {
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


        for(int j=0;j<numVolt;j++) {
            total_charge_normed[j] = GraphIntegral(cc_norm[j],left_edge,right_edge);
        }
        if(total_charge_normed[numVolt-1]<0) {
            for(int j=0;j<numVolt;j++) total_charge_normed[j]=-total_charge_normed[j];
        }
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
    for(int i=0;i<numVolt;i++) sq_volt[i] = sqrt(abs(voltages[i]));

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
    Float_t dv = abs(voltages[1]-voltages[0]);
    depl_fit1->SetRange(sqrt(i_start*dv),sqrt(i_finish*dv));
    depl_fit2->SetRange(sqrt(i_plato*dv),sqrt(voltages[numVolt-1]));

    // Plotting the charge
    TGraph *TotalCg=GraphBuilder(numVolt,sq_volt,total_charge,"Sqrt(abs(Voltage)) [sqrt(V)]","total charge [.arb]","total charge");
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
        dir_depl_normed->cd();

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

        dir_depl->cd();
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
bool ModuleEdgeDepletion::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

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

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}

}


