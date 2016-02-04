/**
 * \file
 * \brief Definition of the TCTReader and TCTWaveform classes. © Particulars d.o.o.
 * \details The code below was originally written by © Particulars d.o.o. and was called TCTAnalyse.
 * If needed one can download the original version from particulars.si.
 *
 */

#include <stdio.h>
#include <iostream>

#ifndef ROOT_TCTReader
#define ROOT_TCTReader
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCT Position Sensitive                                               //
//                                                                      //
// fortmat:                                                             //
// Line 0: Writing type of measurement. 11 = waveform,                  //
// 12 = waveform maximum only, 13 = waveform maximum and integrals      //
// Line 1: Datum                                                        //
// Line 2: Absolute time in s at start                                  //
// Line 2: x0, dx, Nx                                                   //
// Line 3: y0, dy, Ny                                                   //
// Line 4: z0, dz, Nz                                                   //
// Line 5: Nvoltages, Voltages                                          //
// Line 6: t0, dt, Npoints                                              //
// Line a->oo: x y z U tmeas                                            //
// Line b->oo: Npoints meritve                                          //
// Line c->oo: Npoints meritve                                          //
// Line d->oo: Npoints meritve                                          //
// The class TCT scan is used in x                                      //
// The class TCT scan is used in y                                      //
//////////////////////////////////////////////////////////////////////////

#include <TArrayI.h>
#include <TArrayF.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>
class TCTWaveform;

class TCTReader
{
    private:
        void  swoo(char *, char *);
        void  swooip(float *, int);
        int BLE_CODE;
    public:

        FILE *in;
        TArrayI 	 Date;		//Date of the Measurement (or simulation)
        TClonesArray   *histo1;        //->
        TClonesArray   *histo2;        //->
        TClonesArray   *histo3;        //->
        TClonesArray   *histo4;        //->
        Float_t *xyz[9];             // x coordinates coresponing to histo array
        TArrayF ta;                  // ta coordinates coresponing to histo array
        Int_t abstime;
        Int_t Nx,Ny,Nz;
        Int_t RefInd;
        Float_t dx,dy,dz;
        Float_t x0,y0,z0;
        Int_t NU1,NU2;
        TArrayF U1,U2,I1,I2;                    // Array of voltages
        Float_t t0,dt;
        Int_t NP;
        Int_t WFOnOff[4];
        Int_t type;
        Int_t numxyz;
        //Header
        Float_t T;                    // temperature
        Float_t Source;               // type of e-h generation
        Char_t *User;                 // user taking the measurements
        Char_t *Sample;               // Sample name
        Char_t *Comment;              // Comment
        Char_t *FileName;             // The name of the input file

        TCTReader(Char_t *, Float_t=0,Int_t=0);
        ~TCTReader();
        void  ReadWFs(Float_t=0);
        void  ReadWFsBin(Float_t=0);
        Int_t indx(Int_t , Int_t , Int_t, Int_t , Int_t );
        void cords(Int_t, Int_t, Int_t, Int_t, Int_t);
        Float_t GetWidth(TH1F *, Int_t &, Int_t &,  Float_t=0, Float_t=25, Float_t=-1111, Float_t=-1111, Int_t = 50);
        TH1F *GetHA(Int_t , Int_t);
        TH1F *GetHA(Int_t , Int_t, Int_t, Int_t, Int_t=0, Int_t=0);
//        TGraph *GetIV() {TGraph *gr=new TGraph(U1.GetSize(),U1.GetArray(),I1.GetArray()); return gr;}

        TH2F *Draw(Int_t=0, Int_t=0, Int_t=0, Int_t=0, Int_t=0, Int_t=0, Float_t=0, Float_t=25);
        //TH2F *Draw(Int_t=0, Int_t=0, Int_t=0, Int_t=0, Float_t=0, Float_t=25);
        void DrawList(Int_t, Int_t *);
        void DrawList(Int_t, Int_t*, Int_t *);
        void PrintInfo();
        void CorrectBaseLine(Float_t=0.0);

        TCTWaveform *Projection(int ch, int dir,int x,int y,int z, int nu1, int nu2, int num);
        TCTWaveform *Projection(int , int *);

        //ClassDef(TCTReader,1);
};

#endif // ROOT_TCTReader

/////////////////////////////////////////////////////

#ifndef ROOT_TCTWaveform
#define ROOT_TCTWaveform

#include <TArrayI.h>
#include <TArrayF.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>
#include <TPaveText.h>
#include <TGraph.h>

class TCTWaveform {
    private:

        Int_t          Multiple;      //Number of Voltages
        Float_t 	 Frequency;	//Frequency of trigger
        TClonesArray   *histo;        //->
        TArrayF        Voltages;      //Array of Voltages
        TArrayF	 Current; 	//Array of Currents
        TArrayF	 Temperature;	//Array of Temperatures
        TArrayF 	 Date;		//Date of the Measurement (or simulation)
        TArrayF        Frequencies;   //Array of Frequnecies
        Float_t        AnnealTime;

    public:
        TPaveText *pt;                //!

        //stuff added to patch position sensitive TCT
        bool DrawMode;
        Char_t suffix[10];
        Char_t prefix[10];

        TCTWaveform(Int_t=60);
        virtual ~TCTWaveform();

        void AddHisto(Int_t,Float_t,TH1F *); //voltage index, voltage, histogram
        void AddHisto(Float_t,TH1F *); //voltage, histogram
        //Voltage Current Section
        inline void SetAnnealTime(Float_t x) {AnnealTime=x;}
        inline Float_t GetAnnealTime() {return AnnealTime;}
        void SetVoltages(Float_t *voltages) {Voltages.Adopt(Multiple,voltages);};
        void SetVoltage(Int_t index, Float_t voltage) {Voltages[index]=voltage;};
        Int_t int2ascii(Char_t v[],Float_t,Int_t=1);
        void SetHistoTime(int ind,Float_t start);
        void SetHistoTime(int num,Float_t *start);
        Float_t GetVoltage(Int_t index) {return(Voltages[index]);}
        Int_t GetIndex(Float_t Volt) {for(Int_t i=0;i<Multiple;i++) if(Voltages[i]==Volt) return(i); printf("No such measurement!\n"); return(-1);}
        Float_t GetCurrent(Int_t index) {return(Current[index]);}
        Float_t GetCurrent(Float_t Volt) {for(Int_t i=0;i<Multiple;i++) if(Voltages[i]==Volt) return(Current[i]); printf("No such measurement!\n"); return -1;}
        void SetVoltages(Float_t svol,Float_t step) {for(Int_t z=0;z<Multiple;z++) Voltages[z]=svol+step*(Float_t)z;};
        void SetVoltages(Int_t num,Float_t svol,Float_t step) {for(Int_t z=num;z<Multiple;z++) Voltages[z]=svol+step*(Float_t)(z-num);};
        void SetTemperatures(Float_t T)  {for(Int_t z=0;z<Multiple;z++) Temperature[z]=T;};
        void PrintVoltages() {printf("Defined Voltages::\n"); for(Int_t i=0;i<Multiple;i++) printf("(%d,%f)\n",i,Voltages[i]);}

        //Draw Section

        void Draw(Float_t Volt,Option_t *opt,Float_t low=-1111,Float_t high=-1111) {if(GetIndex(Volt)==-1) return; else Draw(GetIndex(Volt),opt,low,high);}; // Voltage, Graphics Option, low bin , high bin
        void DrawTest(int);

        void Draw(Int_t,Option_t *,Float_t=-1111,Float_t=-1111);   // Voltage Index, Graphics Option, low bin , high bin
        void DrawMulti(Float_t=-1111.,Float_t=-1111.,Int_t=-1,Int_t=-1,Int_t =1,Int_t=0,Int_t=-1111); // low bin, high bin, Start Voltage Index, End Voltage Index, Int_t step, Int_t model

        void Legend(TH1 *,Int_t,Int_t,Int_t=1,Int_t =0);

        //
        //Histogram Section
        void GetHistogram(Int_t number,TH1F *his) {((TH1F *)histo->At(number))->Copy(*his);}; //voltage index, histogram to copy original to
        void GetHistogram(Float_t voltage,TH1F *his) {for(Int_t i=0;i<Multiple;i++) if(Voltages[i]==voltage) ((TH1F *)histo->At(i))->Copy(*his);}; //voltage, histogram
        TH1F *GetHA(Float_t voltage)  {Int_t index=0; for(Int_t i=0;i<Multiple;i++) if(Voltages[i]==voltage) index=i; return((TH1F *)histo->At(index));};
        inline TH1F *operator()(Float_t voltage) {return(GetHA(voltage));};
        inline TH1F *operator()(Int_t index) {return((TH1F *)histo->At(index));};
        //TH1F *GetHA(Int_t index) {return((TH1F *)histo->At(index));};
        //Integral Section
        void GetIntegral(Float_t*,Float_t=1,Float_t=-1111,Float_t =-1111); //  in - array of integrals, scale , min time , high time
        void GetIntegral(Float_t*,Float_t*,Float_t=1,Float_t=-1111,Float_t =-1111); // in - array of voltages, in - array of integrals, scale , min time , high time
        Float_t Integral(Float_t,Float_t=-1111,Float_t =-1111); // voltage, mint, maxt of the integration
        Float_t Integral(Int_t,Float_t=-1111,Float_t =-1111); // index of the voltage, mint, maxt of the integration
        TGraph *CCE(Float_t=-1111,Float_t=-1111,Int_t=0,Int_t=1); //Draw CCE: start time, end time, model(sqrt,lin), option (time , T)

        // General section
        Int_t GetEntries() {return(Multiple);};
        void NormArray(Int_t,Float_t *);
        void Info();
        Float_t GetTime(Float_t=3600);
        Float_t GetT() {Float_t sumT=0; for(Int_t i=0;i<Multiple;i++) sumT+=Temperature[i]; return(sumT/Multiple); };
        Float_t GetT(Float_t Volt) {return (Temperature[GetIndex(Volt)]);};
        Float_t GetT(Int_t index) {return (Temperature[index]);};
};

#endif // ROOT_TCTWaveform
