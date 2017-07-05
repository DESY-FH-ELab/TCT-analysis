/**
 * \file
 * \brief Implementation of TCT::ModuleTopFocus methods.
 */

#include "modules/ModuleTopFocus.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPaveStats.h"


namespace TCT {

/// Analyse data (should be reimplemented by developer)
bool ModuleTopFocus::Analysis() {

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
    delete FWHMg;

    // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance

    TGraph *MINQg=GraphBuilder(numO,optical_axis_co,minQ,"optical distance [#mum]","min Q[rel. to max]","Minimum Charge");
    ff_pol3->SetParameter(0,MINQg->GetMean());
    ff_pol3->SetParameter(2,MINQg->GetMinimum());
    MINQg->Fit("my_pol3","q");
    gStyle->SetOptFit(1);
    st = (TPaveStats*)MINQg->FindObject("stats");
    st->SetFitFormat(".5g");
    MINQg->Write("MinCharge");
    delete MINQg;

    // Find the missalignment between z and optical axis
    for(int j=0;j<numO;j++) abs_pos[j]=pos[j]+Sc0;

    TGraph *POSg=GraphBuilder(numO,optical_axis_co,abs_pos,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
    POSg->Fit("pol1","q");
    gStyle->SetOptFit(1);
    POSg->Write("Missalignment");
    delete POSg;

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
        delete FWHMg1;

        // Plotting the minimum charge of the fitted Erfs a.f.o. optical distance with normed data

        TGraph *MINQg1=GraphBuilder(numO,optical_axis_co,minQ_normed,"optical distance [#mum]","min Q[rel. to max]","Minimum Charge");
        ff_pol4->SetParameter(0,MINQg1->GetMean());
        ff_pol4->SetParameter(2,MINQg1->GetMinimum());
        MINQg1->Fit("my_pol4","q");
        gStyle->SetOptFit(1);
        st = (TPaveStats*)MINQg1->FindObject("stats");
        st->SetFitFormat(".5g");
        MINQg1->Write("MinCharge_Normed");
        delete MINQg1;

        // Find the missalignment between z and optical axis
        for(int j=0;j<numO;j++) abs_pos[j]=pos_normed[j]+Sc0;

        TGraph *POSg1=GraphBuilder(numO,optical_axis_co,abs_pos,"optical distance [#mum]","position of the edge [#mum]","Missalignment");
        POSg1->Fit("pol1","q");
        gStyle->SetOptFit(1);
        POSg1->Write("Missalignment_Normed");
        delete POSg1;

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

    delete cc;
    if(config->CH_PhDiode()) delete cc_norm;
    delete ph_charge;
    return true;
}

/// Check before analysis (should be reimplemented by developer)
bool ModuleTopFocus::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

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

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}

}


