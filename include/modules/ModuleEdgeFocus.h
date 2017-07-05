/**
 * \file
 * \brief Definition of the TCT::ModuleEdgeFocus class.
 */

#ifndef __MODULEEdgeFocus_H__
#define __MODULEEdgeFocus_H__ 1

#include "TCTModule.h"

namespace TCT {

  class ModuleEdgeFocus: public TCTModule {

   public:
        ModuleEdgeFocus(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title) {}
        bool CheckModuleData();
        bool Analysis();

    };
}
#endif 
