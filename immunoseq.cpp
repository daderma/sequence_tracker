#include "immunoseq.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <map>


namespace immunoseq
{


std::string get_patient_id(std::string const& sample_id)
{
	std::string result;
	for(auto const& c: sample_id)
	{
		if(std::isdigit(c))
		{
			result += c;
		}
		else
		{
			return result;
		}
	}

	BOOST_THROW_EXCEPTION(load_exception() << load_description_type("unable to extract patient id from sample id: \"" + sample_id + "\""));
}


void load(boost::filesystem::path const& path, patients::patients_type& patients)
{
	char const* sample_name_field("sample_name");
	char const* release_date_field("release_date");
	char const* rearrangement_field("rearrangement");
	char const* rearrangement_type_field("rearrangement_type");
	char const* reads_field("reads");
	char const* v_family_field("v_family");
	char const* d_family_field("d_family");
	char const* j_family_field("j_family");

	boost::filesystem::ifstream stream(path);

	std::map<std::string, std::size_t> header;
	std::string row;
	std::int64_t counter(0);
	while(std::getline(stream, row))
	{
		std::vector<std::string> columns;
		boost::split(columns, row, boost::is_any_of("\t"));
		if(header.empty())
		{
			std::size_t index(0);
			for(auto const& column: columns)
			{
				header.insert(std::make_pair(column, index ++));
			}

			if(header.count(sample_name_field) && header.count(release_date_field)
				&& header.count(rearrangement_field) && header.count(rearrangement_type_field) 
				&& header.count(reads_field) 
				&& header.count(v_family_field) && header.count(d_family_field) && header.count(j_family_field)
				)
			{
				std::cout << "Loading " << path << std::endl;
				continue;
			}
			else
			{
				std::cout << "Ignoring " << path << std::endl;
				break;
			}
		}

		++ counter;
		if(columns.size() != header.size())
		{
			BOOST_THROW_EXCEPTION(load_exception() 
				<< load_description_type("number of columns does not match header")
				<< load_row_type(counter));
		}
		
		auto const rearrangement_type(columns[header[rearrangement_type_field]]);
		if(rearrangement_type != "VDJ")
		{
			/*BOOST_THROW_EXCEPTION(load_exception() 
				<< load_description_type("unsupported rearrangement type: " + columns[header[rearrangement_type_field]])
				<< load_row_type(row));*/
			std::cout << "Warning! Reading rearrangement type \"" << rearrangement_type << "\" at row " << counter << std::endl;
		}

		auto const sample_id(columns[header[sample_name_field]]);
		auto const patient_id(get_patient_id(sample_id));
		auto const timestamp(boost::posix_time::time_from_string(columns[header[release_date_field]]));
		
		auto patient(patients[patient_id]);
		if(!patient)
		{
			patient = std::make_shared<patients::patient_type>();
			patient->id = patient_id;
			patients[patient_id] = patient;
		}

		auto sample(patient->samples[std::make_pair(timestamp, sample_id)]);
		if(!sample)
		{
			sample = std::make_shared<samples::sample_type>();
			sample->id = sample_id;
			sample->timestamp = timestamp;
			patient->samples[std::make_pair(timestamp, sample_id)] = sample;
		}

		auto sequence(std::make_shared<sequences::sequence_type>());
		sequence->rearrangement = columns[header[rearrangement_field]];
		sequence->reads = boost::lexical_cast<std::int64_t>(columns[header[reads_field]]);
		sequence->v_family = columns[header[v_family_field]];
		sequence->d_family = columns[header[d_family_field]];
		sequence->j_family = columns[header[j_family_field]];
		sample->sequences.push_back(sequence);
	}
}


}	// namespace immunoseq