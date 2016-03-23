/**
 * \file
 * \brief Definition of the TCT::ModuleEdgeDepletion class.
 */

#ifndef __MODULEEdgeDepletion_H__
#define __MODULEEdgeDepletion_H__ 1

#include "TCTModule.h"

namespace TCT {

  class ModuleEdgeDepletion: public TCTModule {

   public:
        ModuleEdgeDepletion(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title) {}
        bool CheckModuleData();
        bool Analysis();

    };
}
#endif 
