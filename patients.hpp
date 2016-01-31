#ifndef SEQUENCE_TRACKER__PATIENTS_HPP
#define SEQUENCE_TRACKER__PATIENTS_HPP


#include "samples.hpp"


namespace patients
{


struct patient_type
{
	std::string id;
	samples::samples_type samples;
};


typedef std::shared_ptr<patient_type> patient_ptr_type;
typedef std::map<std::string, patient_ptr_type> patients_type;


}	// namespace patients


#endif	// SEQUENCE_TRACKER__PATIENTS_HPP