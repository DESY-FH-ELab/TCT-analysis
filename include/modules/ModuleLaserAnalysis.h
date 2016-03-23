/**
 * \file
 * \brief Definition of the TCT::ModuleLaserAnalysis class.
 */

#ifndef __MODULELaserAnalysis_H__
#define __MODULELaserAnalysis_H__ 1

#include "TCTModule.h"

namespace TCT {

  class ModuleLaserAnalysis: public TCTModule {

   public:
        ModuleLaserAnalysis(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title) {}
        bool Analysis();
        bool LaserPowerDrop();
        bool BeamSigma();

    };
}
#endif 
