/**
 * \file
 * \brief Definition of the TCT::ModuleTopFocus class.
 */

#ifndef __MODULETOPFocus_H__
#define __MODULETOPFocus_H__ 1

#include "TCTModule.h"

namespace TCT {

  class ModuleTopFocus: public TCTModule {

   public:
        ModuleTopFocus(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title) {}
        bool CheckModuleData();
        bool Analysis();

    };
}
#endif 
