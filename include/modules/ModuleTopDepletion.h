/**
 * \file
 * \brief Definition of the TCT::ModuleTopDepletion class.
 */

#ifndef __MODULETopDepletion_H__
#define __MODULETopDepletion_H__ 1

#include "TCTModule.h"

namespace TCT {

  class ModuleTopDepletion: public TCTModule {

   public:
        ModuleTopDepletion(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title) {}
        bool CheckModuleData();
        bool Analysis();

    };
}
#endif 
