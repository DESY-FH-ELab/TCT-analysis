/**
 * \file
 * \brief Definition of the sample class.
 */

#ifndef __SCANNING_H__
#define __SCANNING_H__ 1

// STD includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TFile.h"
// TCT includes
#include "TCTScan.h"
#include "analysis.h"

namespace TCT {

  class Scanning {

    private :
        int _Nsamples;
    protected:

    public :
        Scanning() {};

        // Default copy constructer should be fine
        Scanning(const Scanning &)               = default;
        Scanning & operator = (const Scanning &) = default;

        // Destructor
        ~Scanning() = default;
        bool ReadTCT(char* filename,analysis* ana, bool HasSubs);
        bool CheckData(PSTCT *stct1, analysis* ana);
        bool DoFocus(TFile* f_rootfile, PSTCT *stct, analysis* ana);
        bool SimulateDoFocus(TFile* f_rootfile, analysis* ana);
        bool CheckFocus(PSTCT *stct1, analysis* ana);


    }; // end of class scanning
}
#endif 
