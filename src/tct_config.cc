/**
 * \file
 * \brief Implementation of TCT::tct_config methods
 */

// STD includes
#include<string>

// TCT includes
#include "tct_config.h"

namespace TCT {

void tct_config::SetParameters(std::map<std::string, std::string> id_val){

    for( auto i : id_val){
        if(i.first == "OutFolder")	_OutFolder = i.second;
        if(i.first == "DataFolder")	_DataFolder = i.second;

        //tct scanning mode parameters
        if(i.first == "CH_Detector")		_CH_Det = atoi((i.second).c_str());
        if(i.first == "CH_Photodiode")	_CH_PhDiode = atoi((i.second).c_str());
        if(i.first == "CH_Trigger")		_CH_Trig = atoi((i.second).c_str());
        if(i.first == "Optical_Axis")		_OptAxis = atoi((i.second).c_str());
        if(i.first == "Focus_Search")		_DO_focus = static_cast<bool>(atoi((i.second).c_str()));
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
        if(i.first == "TopMobility")		_DO_TopMobility = static_cast<bool>(atoi((i.second).c_str()));
        if(i.first == "TopDepletionVoltage")		_DO_TopDepletion = static_cast<bool>(atoi((i.second).c_str()));
        if(i.first == "EdgeDepletionVoltage")		_DO_EdgeDepletion = static_cast<bool>(atoi((i.second).c_str()));
        if(i.first == "EdgeVelocityProfile")		_DO_EdgeVelocity = static_cast<bool>(atoi((i.second).c_str()));
        if(i.first == "EV_Time")          _EV_Time = atof((i.second).c_str());
        if(i.first == "CorrectBias")          _CorrectBias = atof((i.second).c_str());

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

    return;
}

}
