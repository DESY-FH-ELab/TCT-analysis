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

		public :

			param() :
				_MaxAcqs(300),
				_MaxPulses(10),
				_S2nRef(2),
				_S2nCut(3)
			{}

			// no copy constructor needed

			// Dectructor
			~param() = default;

			uint32_t & GetMaxAcqs(){ return _MaxAcqs;}
			const uint32_t & GetMaxAcqs() const { return _MaxAcqs;}

			uint32_t & GetMaxPulses(){ return _MaxPulses;}
			const uint32_t & GetMaxPulses() const { return _MaxPulses;}

			float & GetS2nRef(){ return _S2nRef;}
			const float & GetS2nRef() const { return _S2nRef;}
			
			float & GetS2nCut(){ return _S2nCut;}
			const float & GetS2nCut() const { return _S2nCut;}
			

		};
}

inline std::ostream & operator << (std::ostream & os, const TCT::param & param) {
	return os	<< "\n List of parameters \n"
			<< " MaxAcqs = " 	<< param.GetMaxAcqs()
			<< " MaxPulses = " 	<< param.GetMaxPulses()
			<< " S2nRef = " 	<< param.GetS2nRef()
			<< " S2nCut = " 	<< param.GetS2nCut() << std::endl;

}
#endif
