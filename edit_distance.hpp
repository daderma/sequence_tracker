#ifndef SEQUENCE_TRACKER__EDIT_DISTANCE_HPP
#define SEQUENCE_TRACKER__EDIT_DISTANCE_HPP


#include <cstdint>
#include <string>


namespace edit_distance
{


std::size_t levenshtein(std::string const& left, std::string const & right);


}	// namespace edit_distance


#endif	// SEQUENCE_TRACKER__EDIT_DISTANCE_HPP