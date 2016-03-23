/**
 * \file
 * \brief Definition of the TCT::ModuleEdgeField class.
 */

#ifndef __MODULEEdgeField_H__
#define __MODULEEdgeField_H__ 1

#include "TCTModule.h"

#ifdef USE_GUI
#include "QDoubleSpinBox"
#endif

namespace TCT {

  class ModuleEdgeField: public TCTModule {

   public:
        ModuleEdgeField(tct_config* config1, const char* name, TCT_Type type, const char* title):
            TCTModule(config1, name, type, title), _EV_Time(0.3) {}
        bool CheckModuleData();
        bool Analysis();
#ifdef USE_GUI
        void PrintConfig(std::ofstream &conf_file);
        void AddParameters(QVBoxLayout* layout);
        void FillParameters();
        void ToVariables();
        QDoubleSpinBox* ev_time;
#endif
        void SetEV_Time(float value) { _EV_Time = value; }
        float GetEV_Time() { return _EV_Time; }

  private:
        Double_t ff(Double_t E, Double_t Uuu, Double_t a);
        Double_t BiSectionMethod(Double_t eps, Double_t x1, Double_t x2, Double_t Uuu, Double_t a);
        Double_t Mu(Double_t E, Int_t Type);
        template <typename T> int sgn(T val);
        float _EV_Time;

    };
}
#endif 
