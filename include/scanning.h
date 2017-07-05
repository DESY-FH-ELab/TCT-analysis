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

// TCT includes
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
        bool CreateOutputFile();
        bool Separate_and_Sample();
        bool CheckData();

    }; // end of class scanning
}
#endif 
