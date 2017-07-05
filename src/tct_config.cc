/**
 * \file
 * \brief Implementation of TCT::tct_config methods
 */

// STD includes
#include<string>

// TCT includes
#include "tct_config.h"
#include "modules/ModuleTopFocus.h"
#include "modules/ModuleTopMobility.h"
#include "modules/ModuleTopDepletion.h"
#include "modules/ModuleEdgeFocus.h"
#include "modules/ModuleEdgeDepletion.h"
#include "modules/ModuleEdgeField.h"

namespace TCT {

tct_config::~tct_config() {
    for(int i=0;i<tct_modules.size();i++){
        delete tct_modules[i];
    }
    tct_modules.clear();
}

void tct_config::SetParameters(std::map<std::string, std::string> id_val){

    for( auto i : id_val){
        if(i.first == "OutFolder")	_OutFolder = i.second;
        if(i.first == "DataFolder")	_DataFolder = i.second;

        //tct scanning mode parameters
        if(i.first == "CH_Detector")		_CH_Det = atoi((i.second).c_str());
        if(i.first == "CH_Photodiode")	_CH_PhDiode = atoi((i.second).c_str());
        if(i.first == "CH_Trigger")		_CH_Trig = atoi((i.second).c_str());
        if(i.first == "Optical_Axis")		_OptAxis = atoi((i.second).c_str());
        if(i.first == "Scanning_Axis")    _ScAxis = atoi((i.second).c_str());
        if(i.first == "FWHM")             _FFWHM = atof((i.second).c_str());
        if(i.first == "TimeSensorLow")    _FTlow = atof((i.second).c_str());
        if(i.first == "TimeSensorHigh")   _FThigh = atof((i.second).c_str());
        if(i.first == "TimeDiodeLow")     _FDLow = atof((i.second).c_str());
        if(i.first == "TimeDiodeHigh")    _FDHigh = atof((i.second).c_str());
        if(i.first == "SaveSeparateCharges")		_FSeparateCharges = static_cast<bool>(atoi((i.second).c_str()));
        if(i.first == "SaveSeparateWaveforms")	_FSeparateWaveforms = static_cast<bool>(atoi((i.second).c_str()));
        if(i.first == "Movements_dt")     _Movements_dt = atof((i.second).c_str());
        if(i.first == "TCT_Mode")         _TCT_Mode = atoi((i.second).c_str());
        if(i.first == "Voltage_Source")   _VoltSource = atoi((i.second).c_str());
        if(i.first == "CorrectBias")          _CorrectBias = atof((i.second).c_str());

        //modular system
        if(i.first == "Focus_Search")		{
            RegisterModule(new ModuleTopFocus(this,"Focus_Search",_Top,"Focus Search"),static_cast<bool>(atoi((i.second).c_str())));
            RegisterModule(new ModuleEdgeFocus(this,"Focus_Search",_Edge,"Focus Search"),static_cast<bool>(atoi((i.second).c_str())));
        }
        if(i.first == "TopDepletionVoltage")		RegisterModule(new ModuleTopDepletion(this,"TopDepletionVoltage",_Top,"Depletion Voltage"),static_cast<bool>(atoi((i.second).c_str())));
        if(i.first == "TopMobility")                RegisterModule(new ModuleTopMobility(this,"TopMobility",_Top,"Charge Carriers Mobility"),static_cast<bool>(atoi((i.second).c_str())));
        if(i.first == "EdgeDepletionVoltage")		RegisterModule(new ModuleEdgeDepletion(this,"EdgeDepletionVoltage",_Edge,"Depletion Voltage"),static_cast<bool>(atoi((i.second).c_str())));
        if(i.first == "EdgeVelocityProfile")		RegisterModule(new ModuleEdgeField(this,"EdgeVelocityProfile",_Edge,"Electric Field Profiles"),static_cast<bool>(atoi((i.second).c_str())));
        //end modular system

        //coefficients, constants, parameters
        if(i.first == "Mu0_Electrons")                _mu0_els = atof((i.second).c_str());
        if(i.first == "Mu0_Holes")                    _mu0_holes = atof((i.second).c_str());
        if(i.first == "SaturationVelocity")           _v_sat = atof((i.second).c_str());
        if(i.first == "Amplification")                _ampl = atof((i.second).c_str());
        if(i.first == "LightSplitter")                _light_split = atof((i.second).c_str());
        if(i.first == "ResistanceSensor")             _R_sensor = atof((i.second).c_str());
        if(i.first == "ResistancePhotoDetector")      _R_diode = atof((i.second).c_str());
        if(i.first == "ResponcePhotoDetector")        _RespPhoto = atof((i.second).c_str());
        if(i.first == "EnergyPair")                   _E_pair = atof((i.second).c_str());


    }
    for( auto i : id_val) {
        if(i.first == "EV_Time")                    ((ModuleEdgeField*)GetModule("EdgeVelocityProfile"))->SetEV_Time(atof((i.second).c_str()));
    }

    return;
}

TCTModule* tct_config::GetModule(int index) {
    if(index<tct_modules.size()) return tct_modules.at(index);
    else {
        std::cout<<"!!! TCT Module index out of range.\n";
        return 0;
    }
}

TCTModule* tct_config::GetModule(const char *name) {
    for(int i=0;i<tct_modules.size();i++) {
        if(strcmp(((TCTModule*)tct_modules[i])->GetName(),name)==0) return tct_modules[i];
    }
    std::cout<<"!!! TCT Module not found.\n";
    return 0;

}
void tct_config::RegisterModule(TCTModule* module, bool enabled) {
    module->setEnabled(enabled);
    tct_modules.push_back(module);
}

}
