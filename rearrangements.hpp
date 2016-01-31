#ifndef SEQUENCE_TRACKER__REARRANGEMENTS_HPP
#define SEQUENCE_TRACKER__REARRANGEMENTS_HPP


#include <cstdint>
#include <memory>
#include <map>
#include <string>
#include <ostream>


namespace rearrangements
{


struct rearrangement_type
{
	std::int64_t id;
	std::string nucleotides;
};


typedef std::shared_ptr<rearrangement_type> rearrangement_ptr_type;
typedef std::map<std::string, rearrangement_ptr_type> rearrangements_type;


}	// namespace sequences


#endif	// SEQUENCE_TRACKER__REARRANGEMENTS_HPP