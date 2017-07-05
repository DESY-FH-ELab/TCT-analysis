/**
 * \file
 * \brief Definition of the TCT::TCTModule class.
 */

#ifndef __TCTMODULE_H__
#define __TCTMODULE_H__ 1

// STD includes
#include <iostream>
#include "fstream"
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
// TCT includes
#include "tct_config.h"
#include "TCTReader.h"

#ifdef USE_GUI
#include "QVBoxLayout"
#endif

namespace TCT {

enum TCT_Type {
    _Top = 0,
    _Edge,
    _Bottom
};

  class TCTModule {

    protected:
        bool enabled;
        const char* fName;
        TCT_Type fType;
        const char* fTitle;
        tct_config* config;
        TCTReader* stct;
        TFile* f_rootfile;

    public :
        TCTModule(tct_config* config1, const char* name, TCT_Type type, const char* title);

        // Default copy constructer should be fine
        TCTModule(const TCTModule &)              ;
        TCTModule & operator = (const TCTModule &);

        // getters
        const char* GetName() { return fName; }
        TCT_Type GetType() { return fType; }
        const char* GetTitle() { return fTitle; }
        bool isEnabled() { return enabled; }

        // setters
        void setEnabled(bool value) { enabled = value; }

        // analysis runner
        bool Do(TCTReader *in_stct, TFile *in_rootfile);
        // default analysis
        virtual bool CheckModuleData();
        virtual bool Analysis();

#ifdef USE_GUI
        virtual void PrintConfig(std::ofstream &conf_file) {}
        virtual void AddParameters(QVBoxLayout* layout) {}
        virtual void FillParameters() {}
        virtual void ToVariables() {}
#endif

        // default functions
        void SwitchAxis(Int_t sw, Int_t& nPoints, Float_t& step, Float_t& p0);
        void CalculateCharges(Int_t Channel, Int_t Ax, Int_t numAx,  Int_t scanning, Int_t numS, TGraph **charges, Float_t tstart, Float_t tfinish);
        TGraph** NormedCharge(TGraph** sensor, TGraph** photodiode, Int_t numP);
        void FindEdges(TGraph* gr, Int_t numS, Float_t dx, Double_t &left_edge, Double_t &right_edge);
        void FindEdges(TGraph** gr, Int_t numP, Int_t numS, Float_t dx, Float_t* left_pos, Float_t* left_width, Float_t* right_pos, Float_t* right_width);
        TGraph* GraphBuilder(Int_t N, Float_t *x, Float_t *y, const char *namex, const char *namey, const char *title);
        void GraphBuilder(Int_t N, Float_t *x, Float_t *y, const char *namex, const char *namey, const char *title, const char *write_name);
        TGraph* GraphBuilder(Int_t N, Double_t *x, Double_t *y, const char *namex, const char *namey, const char *title);
        void GraphBuilder(Int_t N, Double_t *x, Double_t *y, const char *namex, const char *namey, const char *title, const char *write_name);
        void GraphSeparate(Int_t N, TGraph **gr, const char *dir_name, const char *namex, const char *namey, const char *title, const char *name_0, Double_t *name_1);
        void GraphSeparate(Int_t N, TGraph **gr, const char *dir_name, const char *namex, const char *namey, const char *title, const char *name_0, Float_t *name_1);
        void MultiGraphWriter(Int_t N, TGraph **gr, const char *namex, const char *namey, const char *title, const char *write_name);
        void SetFitParameters(TF1* ff, Double_t p0, Double_t p1, Double_t p2, Double_t p3);
        Double_t GraphIntegral(TGraph *gr, Double_t x1, Double_t x2);
        Double_t abs(Double_t x);
        void ChargeCorrelationHist(TGraph** sensor, TGraph** photodetector, Int_t numO);

    };
}
#endif 
