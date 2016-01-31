#ifndef SEQUENCE_TRACKER__SAMPLES_HPP
#define SEQUENCE_TRACKER__SAMPLES_HPP


#include "sequences.hpp"
#include <boost/date_time.hpp>
#include <map>


namespace samples
{


struct sample_type
{
	std::string id;
	boost::posix_time::ptime timestamp;
	sequences::sequences_type sequences;
};


typedef std::shared_ptr<sample_type> sample_ptr_type;
typedef std::map<std::pair<boost::posix_time::ptime, std::string>, sample_ptr_type> samples_type;


}	// namespace samples


#endif	// SEQUENCE_TRACKER__SAMPLES_HPP