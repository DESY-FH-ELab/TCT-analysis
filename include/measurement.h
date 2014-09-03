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
      std::string _OutFolder; // Folder, where the root file is writen to !! move to Analysis class
      uint32_t _Nacqs;
      float _Noise_Cut;

    public :
      measurement() = default;

      measurement(std::string folder) : 
        _Name("Meas1"),
        _DataInFolder(folder),
	_Nacqs(5),
	_Noise_Cut(0.004){
      }

      measurement(std::string infolder, std::string outfolder) : 
        _Name("Meas1"),
        _DataInFolder(infolder),
        _OutFolder(outfolder),
	_Nacqs(5),
	_Noise_Cut(0.004){
      }

      // Default copy constructer should be fine
      measurement(const measurement &)               = default;
      measurement & operator = (const measurement &) = default;

      std::string & Folder() {return _DataInFolder;}
      const std::string & Folder() const {return _DataInFolder;}

      std::string & OutFolder() {return _OutFolder;} // !! move to ana class
      const std::string & OutFolder() const {return _OutFolder;}

      std::string & Name() {return _Name;}
      const std::string & Name() const {return _Name;}

      float Noise_Cut() { return _Noise_Cut;}
      void SetNoiseCut(float cut_value) {_Noise_Cut = cut_value;}

      bool AcqsLoader(std::vector<TCT::acquisition_single> *acqs, uint32_t maxAcqs = -1);
      void AcqsAnalyser(TCT::acquisition_single *acq, uint32_t iAcq, TCT::acquisition_avg *acqAvg);
      bool AcqsSelecter(TCT::acquisition_single *acq);
      void AcqsProfileFiller(TCT::acquisition_single *acq, TCT::acquisition_avg *acqAvg);
      void AcqsWriter(TCT::sample* sample, std::vector<TCT::acquisition_single> *acqs, TCT::acquisition_avg *acqAvg);

  };
}

inline std::ostream & operator << (std::ostream & os, const TCT::measurement & meas) {
	return os 	<< "\n This measurement is called \n" << meas.Name()
			<< "   - data folder = " << meas.Folder() << std::endl;
}
#endif
