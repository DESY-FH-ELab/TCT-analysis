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
//#include "laser.h"

//  includes from ROOT libraries

namespace TCT {

  // \brief Measurement class for measurement implementation	

  class measurement {

    private :
      std::string _Name;
      std::string _Folder; // Folder, where the associated data files are
      uint32_t _Nacqs;

    public :
      measurement() = default;

      measurement(std::string folder) : 
        _Name("Meas1"),
        _Folder(folder),
	_Nacqs(5){
      }

      // Default copy constructer should be fine
      measurement(const measurement &)               = default;
      measurement & operator = (const measurement &) = default;

      std::string & Folder() {return _Folder;}
      const std::string & Folder() const {return _Folder;}

      std::string & Name() {return _Name;}
      const std::string & Name() const {return _Name;}

      void AcqsLoader(std::vector<TCT::acquisition_single> acqs, uint32_t maxAcqs = 0);
  };
}

inline std::ostream & operator << (std::ostream & os, const TCT::measurement & meas) {
	return os 	<< "\n This measurement is called \n" << meas.Name()
			<< "   - data folder = " << meas.Folder() << std::endl;
}
#endif
