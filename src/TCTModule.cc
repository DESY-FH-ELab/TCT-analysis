/**
 * \file
 * \brief Implementation of TCT::TCTModule methods.
 */

// STD includes
#include<string>
#include<sstream>
#include<iostream>

// TCT includes
#include "TCTModule.h"
#include "TCTReader.h"

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

namespace TCT {

TCTModule::TCTModule(tct_config* config1, const char* name, TCT_Type type, const char* title):
    config(config1),
    fName(name),
    fType(type),
    fTitle(title)
{
    ;
}

/// DO Method
bool TCTModule::Do(TCTReader *in_stct, TFile *in_rootfile) {

    stct = in_stct;
    f_rootfile = in_rootfile;

    if(!CheckModuleData()) {std::cout<<"Not enough data for "<<GetTitle()<<". Skipping..."<<std::endl; return false;}

    if(Analysis()) {
        std::cout<<GetTitle()<<" finished succesfully!"<<std::endl;
        return true;
    }
    else return false;
}
/// Analyse data (should be reimplemented by developer)
bool TCTModule::Analysis() {
    stct->PrintInfo();
    return true;
}

/// Check before analysis (should be reimplemented by developer)
bool TCTModule::CheckModuleData() {

    std::cout<<"\t- DO "<<GetTitle()<<"(type - "<<GetType()<<")\n";

    //check the data here

    std::cout<<GetTitle()<<" Data Test Passed. Processing..."<<std::endl;
    return true;
}


/// To set number of scanning points, step and x0 for each axis.
void TCTModule::SwitchAxis(Int_t sw, Int_t& nPoints, Float_t& step, Float_t& p0) {
    switch(sw)
    {
    case 0: nPoints=stct->Nx; step=stct->dx; p0=stct->x0;  break; //x-axis
    case 1: nPoints=stct->Ny; step=stct->dy; p0=stct->y0; break; //y-axis
    case 2: nPoints=stct->Nz; step=stct->dz; p0=stct->z0; break; //z-axis
    }
}

/// Calculate Charges for given Waveforms
void TCTModule::CalculateCharges(Int_t Channel, Int_t Ax, Int_t numAx, Int_t scanning, Int_t numS, TGraph **charges, Float_t tstart, Float_t tfinish) {
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
    delete wf;
}

/// Calculate Normed Charges
TGraph** TCTModule::NormedCharge(TGraph** sensor, TGraph** photodiode, Int_t numP) {

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

    delete temp_normed;
    return normed_charge;

}

/// Find position of two edge for a fixed voltage (at least fully depleted) and for one optical distance
void TCTModule::FindEdges(TGraph* gr, Int_t numS, Float_t dx, Double_t& left_edge, Double_t& right_edge) {

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

    delete ff_left;
    delete ff_right;
    delete yy;

}

/// Find position of two edge for a fixed voltage (at least fully depleted) and for all optical distances
void TCTModule::FindEdges(TGraph** gr, Int_t numP, Int_t numS, Float_t dx, Float_t* left_pos, Float_t* left_width, Float_t* right_pos, Float_t* right_width) {

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

/// Build graph and return pointer
TGraph* TCTModule::GraphBuilder(Int_t N, Float_t* x, Float_t* y,const char* namex,const char* namey, const char* title) {
    TGraph *temp=new TGraph(N,x,y);
    temp->SetMarkerStyle(21);
    temp->Draw("AP");
    temp->GetHistogram()->GetXaxis()->SetTitle(namex);
    temp->GetHistogram()->GetYaxis()->SetTitle(namey);
    temp->SetTitle(title);
    return temp;
}

/// Build graph, write to file, and clean memory
void TCTModule::GraphBuilder(Int_t N, Float_t* x, Float_t* y,const char* namex,const char* namey, const char* title, const char* write_name) {
    TGraph *temp=new TGraph(N,x,y);
    temp->SetMarkerStyle(21);
    temp->Draw("AP");
    temp->GetHistogram()->GetXaxis()->SetTitle(namex);
    temp->GetHistogram()->GetYaxis()->SetTitle(namey);
    temp->SetTitle(title);
    temp->Write(write_name);
    delete temp;
    //return temp;
}

/// Build graph and return pointer
TGraph* TCTModule::GraphBuilder(Int_t N, Double_t* x, Double_t* y,const char* namex,const char* namey, const char* title) {
    TGraph *temp=new TGraph(N,x,y);
    temp->SetMarkerStyle(21);
    temp->Draw("AP");
    temp->GetHistogram()->GetXaxis()->SetTitle(namex);
    temp->GetHistogram()->GetYaxis()->SetTitle(namey);
    temp->SetTitle(title);
    return temp;
}

/// Build graph, write to file, and clean memory
void TCTModule::GraphBuilder(Int_t N, Double_t* x, Double_t* y,const char* namex,const char* namey, const char* title, const char* write_name) {
    TGraph *temp=new TGraph(N,x,y);
    temp->SetMarkerStyle(21);
    temp->Draw("AP");
    temp->GetHistogram()->GetXaxis()->SetTitle(namex);
    temp->GetHistogram()->GetYaxis()->SetTitle(namey);
    temp->SetTitle(title);
    temp->Write(write_name);
    delete temp;
    //return temp;
}

/// Write separate waveforms to the file
void TCTModule::GraphSeparate(Int_t N, TGraph **gr, const char *dir_name, const char *namex, const char *namey, const char *title, const char *name_0, Double_t *name_1) {

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

/// Write separate waveforms to the file
void TCTModule::GraphSeparate(Int_t N, TGraph **gr, const char *dir_name, const char *namex, const char *namey, const char *title, const char *name_0, Float_t *name_1) {

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

/// Write multigraph to the file
void TCTModule::MultiGraphWriter(Int_t N, TGraph **gr, const char *namex, const char *namey, const char *title, const char *write_name) {

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

    delete mg;

}

/// Brief setting of the fit parameters.
void TCTModule::SetFitParameters(TF1* ff,Double_t p0,Double_t p1,Double_t p2,Double_t p3) {
    ff->SetParameter(0,p0);
    ff->SetParameter(1,p1);
    ff->SetParameter(2,p2);
    ff->SetParameter(3,p3);
}

/// Integrate graph in the interval.
Double_t TCTModule::GraphIntegral(TGraph *gr, Double_t x1, Double_t x2) {
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

/// Quick abs implementation
Double_t TCTModule::abs(Double_t x) {
    if(x<0) return -x;
    return x;
}

/// Correlation between sensor and photo-diode charge
void TCTModule::ChargeCorrelationHist(TGraph** sensor, TGraph** photodetector, Int_t numO) {

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

}


