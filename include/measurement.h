/**
 * \file
 * \brief Definition of the measurement class.
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
      uint32_t _Nacqs;

    public :

      measurement(std::string infolder) : 
        measurement ()
      {
        SetDataInFolder(infolder);
      }

      measurement() :
        _Name("Meas1"),
	_Nacqs(5)
      {}

      // Default copy constructer should be fine
      measurement(const measurement &)               = default;
      measurement & operator = (const measurement &) = default;

      std::string DataInFolder() {return _DataInFolder;}
      void SetDataInFolder(std::string string) { _DataInFolder = string;}
      const std::string & DataInFolder() const {return _DataInFolder;}

      std::string Name() {return _Name;}
      const std::string & Name() const {return _Name;}

      bool AcqsLoader(std::vector<TCT::acquisition_single> *acqs, uint32_t maxAcqs = -1);
      //void AcqsAnalyser(TCT::acquisition_single *acq, uint32_t iAcq, TCT::acquisition_avg *acqAvg);
      //bool AcqsSelecter(TCT::acquisition_single *acq);
      //void AcqsProfileFiller(TCT::acquisition_single *acq, TCT::acquisition_avg *acqAvg);
      //void AcqsWriter(TCT::sample* sample, std::vector<TCT::acquisition_single> *acqs, TCT::acquisition_avg *acqAvg);

  };
}

inline std::ostream & operator << (std::ostream & os, const TCT::measurement & meas) {
	return os 	<< "\n This measurement is called \n" << meas.Name()
			<< "   - data folder = " << meas.DataInFolder() << std::endl;
}
#endif
