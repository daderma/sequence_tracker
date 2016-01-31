#ifndef SEQUENCE_TRACKER__IMMUNOSEQ_HPP
#define SEQUENCE_TRACKER__IMMUNOSEQ_HPP


#include "patients.hpp"
#include <boost/filesystem.hpp>
#include <boost/exception/all.hpp>


namespace immunoseq
{


struct load_exception : virtual std::exception, virtual boost::exception {};
typedef boost::error_info<struct load_description_, std::string> load_description_type;
typedef boost::error_info<struct load_row_, std::int64_t> load_row_type;


void load(boost::filesystem::path const& path, patients::patients_type& patients);


}	// namespace immunoseq


#endif	// SEQUENCE_TRACKER__IMMUNOSEQ_HPP