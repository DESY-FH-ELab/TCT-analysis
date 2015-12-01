/**
 * \file
 * \brief Definition of the TCT::sample class.
 */

#ifndef __SAMPLE_H__
#define __SAMPLE_H__ 1

//  STD includes
#include <iostream>
#include <string>
#include <map>

//  TCT includes
//#include "pulse.h"
//#include "laser.h"


namespace TCT {

  // \brief Sample class for sample implementation	

  class sample {

    private :
      std::string _SampleID;	// identifier for sample
      std::string _Folder;	// folder where sample card is stored
      double _Thickness;


    public :

      // the empty constructor sets a dummy strong and a dummy thickness
      sample() : 
	_SampleID("def"),
	_Thickness(-1){}


      sample(std::string folder) : 
	_Folder(folder),
	_Thickness(-1){}

      sample(std::map<std::string, std::string> id_val) : 
	sample()
      {
        SetParameters(id_val);        
      }

      // Default copy constructer should be fine
      sample(const sample &)               = default;
      sample & operator = (const sample &) = default;

      // Dectructor
      ~sample() = default;

      double Thickness(){ return _Thickness;}
      void SetThickness(float val) {_Thickness = val;}
      const double & Thickness() const{ return _Thickness;}

      std::string SampleID() {return _SampleID;}
      void SetSampleID(std::string val) {_SampleID = val;}
      const std::string & SampleID() const {return _SampleID;}

      //void SetThickness(double thick) { _Thickness = thick;}
      //void SetSampleID(std::string sampleID) { _SampleID = sampleID;}

      std::string Folder() { return _Folder;}
      void SetFolder(std::string folder) {_Folder = folder;}
      const std::string & Folder() const { return _Folder;}

      void SetParameters(std::map<std::string, std::string> id_val);

  };

}

inline std::ostream & operator << (std::ostream & os, const TCT::sample & sam) {
  return os 	<< "\n This sample is called " << sam.SampleID()
    << "\n   - SampleCardFolder = " << sam.Folder() 
    << "\n   - Sample ID = " << sam.SampleID() 
    << "\n   - thickness = " << sam.Thickness() 
    << std::endl;
}
#endif
