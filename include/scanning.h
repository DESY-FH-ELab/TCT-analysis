/**
 * \file
 * \brief Definition of the TCT::Scanning class.
 */

#ifndef __SCANNING_H__
#define __SCANNING_H__ 1

// STD includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
// TCT includes
//#include "TCTScan.h"
#include "tct_config.h"

#ifdef USE_GUI
#include "gui_consoleoutput.h"
#endif

class TCTReader;

namespace TCT {

  class Scanning {

    private :
        int _Nsamples;
        TFile* f_rootfile;
        TCTReader* stct;
        tct_config* config;

    protected:

    public :
        Scanning() {};

        // Default copy constructer should be fine
        Scanning(const Scanning &)              ;
        Scanning & operator = (const Scanning &);

        // Destructor
        //~Scanning();
#ifndef USE_GUI
        bool ReadTCT(char* filename, tct_config* config1);
#else
        bool ReadTCT(char* filename, tct_config* config1, Ui::ConsoleOutput *progress);
#endif
        bool DoTopFocus();
        bool DoTopDepletion();
        bool DoTopMobility();
        bool DoEdgeFocus();
        bool DoEdgeDepletion();
        bool DoEdgeVelocity();
        bool LaserPowerDrop();
        bool BeamSigma();
        bool CheckData();
        bool CheckFocus();
        bool CheckTopDepletion();
        bool CheckTopMobility();
        bool CheckEdgeDepletion();
        bool CheckEdgeVelocity();
        void SwitchAxis(Int_t sw, Int_t& nPoints, Float_t& step, Float_t& p0);
        void CalculateCharges(Int_t Channel, Int_t Ax, Int_t numAx,  Int_t scanning, Int_t numS, TGraph **charges, Float_t tstart, Float_t tfinish);
        TGraph** NormedCharge(TGraph** sensor, TGraph** photodiode, Int_t numP);
        void ChargeCorrelationHist(TGraph** sensor, TGraph** photodetector, Int_t numO);
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
        Double_t ff(Double_t E, Double_t Uuu, Double_t a);
        Double_t BiSectionMethod(Double_t eps, Double_t x1, Double_t x2, Double_t Uuu, Double_t a);
        Double_t Mu(Double_t E, Int_t Type);
        template <typename T> int sgn(T val);
        Double_t abs(Double_t x);

    }; // end of class scanning
}
#endif 
