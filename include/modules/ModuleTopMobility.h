/**
 * \file
 * \brief Definition of the TCT::ModuleTopMobility class.
 */

#ifndef __MODULETOPMOBILITY_H__
#define __MODULETOPMOBILITY_H__ 1

#include "TCTModule.h"

namespace TCT {

  class ModuleTopMobility: public TCTModule {

   public:
        ModuleTopMobility(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title) {}
        bool CheckModuleData();
        bool Analysis();

    };
}
#endif 
