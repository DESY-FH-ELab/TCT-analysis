// Standard includes
#include <iostream>
#include <vector>
#include <map>

namespace TCT {

#ifndef TCT_CONFIG_H
#define TCT_CONFIG_H

class tct_config
{

private :

    uint32_t _TCT_Mode;
    float _SampleThickness;

    std::string _OutFolder; // Folder, where the root file is writen to
    std::string _OutSample_ID;
    std::string _DataFolder; //

    //scanning parameters
    uint32_t _CH_Det;
    uint32_t _CH_PhDiode;
    uint32_t _CH_Trig;
    uint32_t _OptAxis;
    uint32_t _VoltSource;
    bool _DO_focus;
    bool _DO_EdgeDepletion;
    bool _DO_EdgeVelocity;
    uint32_t _ScAxis;
    float _FFWHM;
    float _FTlow;
    float _FThigh;
    float _FDLow;
    float _FDHigh;
    bool _FSeparateCharges;
    bool _FSeparateWaveforms;
    float _Movements_dt;
    float _EV_Time;

    //coefficients
    float _mu0_els;
    float _mu0_holes;
    float _v_sat;
    float _ampl;
    float _light_split;
    float _R_sensor;
    float _R_diode;
    float _RespPhoto;
    float _E_pair;

public:

    tct_config() :
        _OutFolder("../results/"),
        _DataFolder("def"),
        _CH_Det(0),
        _CH_PhDiode(0),
        _CH_Trig(0),
        _OptAxis(3),
        _DO_focus(false),
        _DO_EdgeDepletion(false),
        _DO_EdgeVelocity(false),
        _Movements_dt(0),
        _EV_Time(0.3),
        _FFWHM(10.),
        _TCT_Mode(0)
    {
        //std::cout << "\n   *** No parameter map passes, using default cut values! ***" << std::endl;
    }


    tct_config(std::map<std::string, std::string> id_val) :
        tct_config()
    {
        SetParameters(id_val);
    }


    uint32_t TCT_Mode() { return _TCT_Mode;}
    void SetTCT_Mode(uint32_t val) { _TCT_Mode = val;}
    const uint32_t & TCT_Mode() const { return _TCT_Mode;}

    //begin scanning section
    uint32_t CH_Det() { return _CH_Det;}
    void SetCH_Det(uint32_t val) { _CH_Det = val;}
    const uint32_t & CH_Det() const { return _CH_Det;}

    uint32_t CH_PhDiode() { return _CH_PhDiode;}
    void SetCH_PhDiode(uint32_t val) { _CH_PhDiode = val;}
    const uint32_t & CH_PhDiode() const { return _CH_PhDiode;}

    uint32_t CH_Trig() { return _CH_Trig;}
    void SetCH_Trig(uint32_t val) { _CH_Trig = val;}
    const uint32_t & CH_Trig() const { return _CH_Trig;}

    uint32_t OptAxis() { return _OptAxis;}
    void SetOptAxis(uint32_t val) { _OptAxis = val;}
    const uint32_t & OptAxis() const { return _OptAxis;}

    uint32_t VoltSource() { return _VoltSource;}
    void SetVoltSource(uint32_t val) { _VoltSource = val;}
    const uint32_t & VoltSource() const { return _VoltSource;}

    bool DO_focus() { return _DO_focus;}
    void SetDO_focus(bool val) { _DO_focus = val;}
    const bool & DO_focus() const { return _DO_focus;}

    bool DO_EdgeDepletion() { return _DO_EdgeDepletion;}
    void SetDO_EdgeDepletion(bool val) { _DO_EdgeDepletion = val;}
    const bool & DO_EdgeDepletion() const { return _DO_EdgeDepletion;}

    bool DO_EdgeVelocity() { return _DO_EdgeVelocity;}
    void SetDO_EdgeVelocity(bool val) { _DO_EdgeVelocity = val;}
    const bool & DO_EdgeVelocity() const { return _DO_EdgeVelocity;}

    uint32_t ScAxis() { return _ScAxis;}
    void SetScAxis(uint32_t val) { _ScAxis = val;}
    const uint32_t & ScAxis() const { return _ScAxis;}

    float FFWHM() { return _FFWHM;}
    void SetFFWHM(float val) { _FFWHM = val;}
    const float & FFWHM() const { return _FFWHM;}

    float FTlow() { return _FTlow;}
    void SetFTlow(float val) { _FTlow = val;}
    const float & FTlow() const { return _FTlow;}

    float FThigh() { return _FThigh;}
    void SetFThigh(float val) { _FThigh = val;}
    const float & FThigh() const { return _FThigh;}

    float FDLow() { return _FDLow;}
    void SetFDLow(float val) { _FDLow = val;}
    const float & FDLow() const { return _FDLow;}

    float FDHigh() { return _FDHigh;}
    void SetFDHigh(float val) { _FDHigh = val;}
    const float & FDHigh() const { return _FDHigh;}

    bool FSeparateCharges() { return _FSeparateCharges;}
    void SetFSeparateCharges(bool val) { _FSeparateCharges = val;}
    const bool & FSeparateCharges() const { return _FSeparateCharges;}

    bool FSeparateWaveforms() { return _FSeparateWaveforms;}
    void SetFSeparateWaveforms(bool val) { _FSeparateWaveforms = val;}
    const bool & FSeparateWaveforms() const { return _FSeparateWaveforms;}

    float Movements_dt() { return _Movements_dt;}
    void SetMovements_dt(float val) { _Movements_dt = val;}
    const float & Movements_dt() const { return _Movements_dt;}

    float EV_Time() { return _EV_Time;}
    void SetEV_Time(float val) { _EV_Time = val;}
    const float & EV_Time() const { return _EV_Time;}

    //end scanning section

    //begin coeffiecients

    float mu0_els() { return _mu0_els; }
    void Setmu0_els(float val) { _mu0_els = val;}
    float mu0_holes() { return _mu0_holes; }
    void Setmu0_holes(float val) { _mu0_holes = val;}
    float v_sat() { return _v_sat; }
    void Setv_sat(float val) { _v_sat = val;}
    float ampl() { return _ampl; }
    void Setampl(float val) { _ampl = val;}
    float light_split() { return _light_split; }
    void Setlight_split(float val) { _light_split = val;}
    float R_sensor() { return  _R_sensor; }
    void SetR_sensor(float val) { _R_sensor = val;}
    float R_diode() { return  _R_diode; }
    void SetR_diode(float val) { _R_diode = val;}
    float RespPhoto() { return  _RespPhoto; }
    void SetRespPhoto(float val) { _RespPhoto = val;}
    float E_pair() { return  _E_pair; }
    void SetE_pair(float val) { _E_pair = val;}

    //end coefficients

    std::string DataFolder () { return _DataFolder;}
    void SetDataFolder(std::string val) { _DataFolder = val;}
    const std::string & DataFolder() const { return _DataFolder;}

    void SetParameters(std::map<std::string, std::string> id_val);

    std::string OutFolder() {return _OutFolder;}
    void SetOutFolder(std::string string) { _OutFolder = string;}
    const std::string & OutFolder() const {return _OutFolder;}

    std::string OutSample_ID() {return _OutSample_ID;}
    void SetOutSample_ID(std::string string) { _OutSample_ID = string;}
    const std::string & OutSample_ID() const {return _OutSample_ID;}

    float SampleThickness() { return _SampleThickness;}
    void SetSampleThickness(float val) { _SampleThickness = val;}
    const float & SampleThickness() const { return _SampleThickness;}

};


#endif // TCT_CONFIG_H

#ifndef MODE_SELECTOR_H
#define MODE_SELECTOR_H

class mode_selector
{

private :

    uint32_t _Mode;
    std::string _SampleCard;

public:

    mode_selector() :
        _Mode(0),
        _SampleCard("def")
    {
        //std::cout << "\n   *** No parameter map passes, using default cut values! ***" << std::endl;
    }


    mode_selector(std::map<std::string, std::string> id_val) :
        mode_selector()
    {
        SetParameters(id_val);
    }


    uint32_t Mode() { return _Mode;}
    void SetMode(uint32_t val) { _Mode = val;}
    const uint32_t & Mode() const { return _Mode;}

    std::string SampleCard () { return _SampleCard;}
    void SetSampleCard(std::string val) { _SampleCard = val;}
    const std::string & SampleCard() const { return _SampleCard;}

    void SetParameters(std::map<std::string, std::string> id_val) {
        for( auto i : id_val){
          if(i.first == "Mode") _Mode = atoi((i.second).c_str());
          if(i.first == "SampleCard")	_SampleCard = i.second;
        }
        return;
    }

};


#endif // MODE_SELECTOR_H

}
