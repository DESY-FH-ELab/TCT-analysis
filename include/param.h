/**
 * \file
 * \brief Definition of the parameter class.
 */

#ifndef __PARAM_H__
#define __PARAM_H__ 1

// Standard includes
#include<stdint.h>

namespace TCT {

		// \brief Parameter class for setting of various parameters in the code

		class param {

		private :
			uint32_t _MaxAcqs;
			uint32_t _MaxPulses;

			float _S2nRef;
			float _S2nCut;

			float _Polarity;

		public :

			param() :
				_MaxAcqs(4),
				_MaxPulses(10),
				_S2nRef(2),
				_S2nCut(3)
			{}

			// no copy constructor needed

			// Dectructor
			~param() = default;

			uint32_t & MaxAcqs(){ return _MaxAcqs;}
			const uint32_t & MaxAcqs() const { return _MaxAcqs;}

			uint32_t & MaxPulses(){ return _MaxPulses;}
			const uint32_t & MaxPulses() const { return _MaxPulses;}

			float & S2nRef(){ return _S2nRef;}
			const float & S2nRef() const { return _S2nRef;}
			
			float & S2nCut(){ return _S2nCut;}
			const float & S2nCut() const { return _S2nCut;}

			float & Polarity(){ return _Polarity;}
			const float & Polarity() const { return _Polarity;}		

			//float & (){ return _;}
			//const float & () const { return _;}		

		};
}

inline std::ostream & operator << (std::ostream & os, const TCT::param & param) {
	return os	<< "\n List of parameters \n"
			<< " MaxAcqs = " 	<< param.MaxAcqs()
			<< " MaxPulses = " 	<< param.MaxPulses()
			<< " S2nRef = " 	<< param.S2nRef()
			<< " S2nCut = " 	<< param.S2nCut() << std::endl;

}
#endif
