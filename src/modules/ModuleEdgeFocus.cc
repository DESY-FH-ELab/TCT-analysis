/**
 * \file
 * \brief Implementation of TCT::ModuleEdgeFocus methods.
 */

#include "modules/ModuleEdgeFocus.h"

#include "TStyle.h"
#include "TPaveStats.h"

namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleEdgeFocus::Analysis() {

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

    TGraph *FWHMg_Left = GraphBuilder(numO,optical_axis_co,width_left,"optical distance [#mum]","FWHM [#mum]","Gaussian Beam Profile");
    ff_pol1->SetParameter(0,FWHMg_Left->GetMean());
    ff_pol1->SetParameter(2,FWHMg_Left->GetMinimum());
    FWHMg_Left->Fit("my_pol1","q");
    gStyle->SetOptFit(1);
    st = (TPaveStats*)FWHMg_Left->FindObject("stats");
    st->SetFitFormat(".5g");
    FWHMg_Left->Write("FWHM_Left");
    delete FWHMg_Left;

    TGraph *FWHMg_Right = GraphBuilder(numO,optical_axis_co,width_right,"optical distance [#mum]","FWHM [#mum]","Gaussian Beam Profile");
    ff_pol1->SetParameter(0,FWHMg_Right->GetMean());
    ff_pol1->SetParameter(2,FWHMg_Right->GetMinimum());
    FWHMg_Right->Fit("my_pol1","q");
    gStyle->SetOptFit(1);
    st = (TPaveStats*)FWHMg_Right->FindObject("stats");
    st->SetFitFormat(".5g");
    FWHMg_Right->Write("FWHM_Right");
    delete FWHMg_Right;

    // Find the missalignment between z and optical axis
    for(int j=0;j<numO;j++) {
        abs_pos_left[j]=pos_left[j]+Sc0;
        abs_pos_right[j]=pos_right[j]+Sc0;
    }

    TGraph *POSg_Left = GraphBuilder(numO,optical_axis_co,abs_pos_left,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
    POSg_Left->Fit("pol1","q");
    gStyle->SetOptFit(1);
    POSg_Left->Write("Missalignment_Left");
    delete POSg_Left;

    TGraph *POSg_Right = GraphBuilder(numO,optical_axis_co,abs_pos_right,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
    POSg_Right->Fit("pol1","q");
    gStyle->SetOptFit(1);
    POSg_Right->Write("Missalignment_Right");
    delete POSg_Right;

    // Plotting the sensor thickness
    GraphBuilder(numO,optical_axis_co,sensor_thick,"optical distance [#mum]","sensor thickness [#mum]","Sensor Thickness","SensorThickness");


    if(config->CH_PhDiode()) {
        dir_fsearch_normed->cd();

        //draw the gaussian beam profile with normed charge
        TGraph *FWHMg1_Left = GraphBuilder(numO,optical_axis_co,width_left_normed,"optical distance [#mum]","FWHM [#mum]","Gaussian Beam Profile");
        ff_pol2->SetParameter(0,FWHMg1_Left->GetMean());
        ff_pol2->SetParameter(2,FWHMg1_Left->GetMinimum());
        FWHMg1_Left->Fit("my_pol2","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg1_Left->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg1_Left->Write("FWHM_Left_Normed");
        delete FWHMg1_Left;

        //draw the gaussian beam profile with normed charge
        TGraph *FWHMg1_Right = GraphBuilder(numO,optical_axis_co,width_right_normed,"optical distance [#mum]","FWHM [#mum]","Gaussian Beam Profile");
        ff_pol2->SetParameter(0,FWHMg1_Right->GetMean());
        ff_pol2->SetParameter(2,FWHMg1_Right->GetMinimum());
        FWHMg1_Right->Fit("my_pol2","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)FWHMg1_Right->FindObject("stats");
        st->SetFitFormat(".5g");
        FWHMg1_Right->Write("FWHM_Right_Normed");
        delete FWHMg1_Right;

        // Find the missalignment between z and optical axis
        for(int j=0;j<numO;j++) {
            abs_pos_left[j]=pos_left_normed[j]+Sc0;
            abs_pos_right[j]=pos_right_normed[j]+Sc0;
        }

        TGraph *POSg1_Left = GraphBuilder(numO,optical_axis_co,abs_pos_left,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
        POSg1_Left->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg1_Left->Write("Missalignment_Left_Normed");
        delete POSg1_Left;

        TGraph *POSg1_Right = GraphBuilder(numO,optical_axis_co,abs_pos_right,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
        POSg1_Right->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg1_Right->Write("Missalignment_Right_Normed");
        delete POSg1_Right;

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

    delete cc;
    if(config->CH_PhDiode()) delete cc_norm;
    delete ph_charge;

    return true;
}

/// Check before analysis (should be reimplemented by developer)
bool ModuleEdgeFocus::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

    Int_t NOpt;
    Int_t NSc;
    std::cout<<"\t\t-- Optical Axis: "<< config->OptAxis() <<std::endl;
    switch(config->OptAxis()) {
    case 1: NOpt = stct->Nx; break;
    case 2: NOpt = stct->Ny; break;
    case 3: NOpt = stct->Nz; break;
    }
    if(NOpt>=3) std::cout<<"\t\t\tOptical axis scan contains "<<NOpt<<" points. OK"<<std::endl;
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

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}

}


