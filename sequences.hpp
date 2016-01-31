#ifndef SEQUENCE_TRACKER__SEQUENCES_HPP
#define SEQUENCE_TRACKER__SEQUENCES_HPP


#include "rearrangements.hpp"
#include <vector>


namespace sequences
{


struct sequence_type
{
	std::int64_t id;
	rearrangements::rearrangement_ptr_type rearrangement;
	std::size_t reads;
	std::string v_family;
	std::string d_family;
	std::string j_family;
};


typedef std::shared_ptr<sequence_type> sequence_ptr_type;
typedef std::vector<sequence_ptr_type> sequences_type;


std::ostream& operator<<(std::ostream& stream, sequence_ptr_type const& sequence);


}	// namespace sequences


#endif	// SEQUENCE_TRACKER__SEQUENCES_HPP