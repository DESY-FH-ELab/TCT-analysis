/**
 * \file
 * \brief Definition of the TCT::measurement class.
 */

#ifndef __MEASUREMENT_H__
#define __MEASUREMENT_H__ 1

//  includes from standard libraries
#include <iostream>
#include<string>

//  includes from TCT classes
#include "acquisition.h"
#include "sample.h"
//#include "laser.h"

//  includes from ROOT libraries

namespace TCT {

  // \brief Measurement class for measurement implementation	

  class measurement {

    private :
      std::string _Name;
      std::string _DataInFolder; // Folder, where the associated data files are

      bool _IsManyPulseStructure; //
      bool _IsMIP; //

      float _Temp;	// 
      float _BiasVolt;	// 

      double _X, _Y; // in um. ?? do we need local and global coordinates?

    public :

      measurement(std::string infolder) : 
	measurement (){
	  SetDataInFolder(infolder);
	}

      measurement() :
	_Name("Meas1"),
	_IsMIP(false),
	_Temp(295.),
	_BiasVolt(500.){
	}

      // Default copy constructer should be fine
      measurement(const measurement &)              ;
      measurement & operator = (const measurement &);

      std::string DataInFolder() {return _DataInFolder;}
      void SetDataInFolder(std::string string) { _DataInFolder = string;}
      const std::string & DataInFolder() const {return _DataInFolder;}

      std::string Name() {return _Name;}
      const std::string & Name() const {return _Name;}

      void SetBiasVolt(float volt) { _BiasVolt = volt;} 
      float BiasVolt() {return _BiasVolt;} 
      const float & BiasVolt() const { return _BiasVolt;} 

      void SetTemp(float temp) { _Temp = temp;} 
      float Temp() {return _Temp;}  
      const float & Temp() const { return _Temp;} 

      bool IsMIP() { return _IsMIP;} 
      void SetIsMIP(bool Is) { _IsMIP = Is;} 

      bool AcqsLoader(std::vector<TCT::acquisition_single> *acqs, uint32_t maxAcqs = -1, bool LeCroyRAW = false);

  };
}

inline std::ostream & operator << (std::ostream & os, const TCT::measurement & meas) {
  return os 	<< "\n This measurement is called \n" << meas.Name()
    << "   - data folder = " << meas.DataInFolder() 
    << std::endl;
}
#endif
