#ifndef SEQUENCE_TRACKER__SEQUENCES_HPP
#define SEQUENCE_TRACKER__SEQUENCES_HPP


#include <cstdint>
#include <memory>
#include <vector>
#include <string>


namespace sequences
{


struct sequence_type
{
	std::string rearrangement;
	std::size_t reads;
	std::string v_family;
	std::string d_family;
	std::string j_family;
};


typedef std::shared_ptr<sequence_type> sequence_ptr_type;
typedef std::vector<sequence_ptr_type> sequences_type;


}	// namespace sequences


#endif	// SEQUENCE_TRACKER__SEQUENCES_HPP