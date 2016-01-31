#include "immunoseq.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <map>


namespace immunoseq
{


void load(boost::filesystem::path const& path, samples::samples_type& samples)
{
	char const* sample_name_field("sample_name");
	char const* sample_tags_field("sample_tags");
	char const* rearrangement_field("rearrangement");
	char const* rearrangement_type_field("rearrangement_type");
	char const* reads_field("reads");
	char const* v_family_field("v_family");
	char const* d_family_field("d_family");
	char const* j_family_field("j_family");

	boost::filesystem::ifstream stream(path);

	rearrangements::rearrangements_type rearrangements;
	std::int64_t sequence_id(0);
		
	std::map<std::string, std::size_t> header;
	std::string row;
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

			if(header.count(sample_name_field) && header.count(sample_tags_field) 
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

		if(columns.size() != header.size())
		{
			BOOST_THROW_EXCEPTION(load_exception() 
				<< load_description_type("number of columns does not match header")
				<< load_row_type(sequence_id));
		}
		
		auto const rearrangement_type(columns[header[rearrangement_type_field]]);
		if(rearrangement_type != "VDJ")
		{
			/*BOOST_THROW_EXCEPTION(load_exception() 
				<< load_description_type("unsupported rearrangement type: " + columns[header[rearrangement_type_field]])
				<< load_row_type(row));*/
		}

		auto& sample(samples[columns[header[sample_name_field]]]);
		if(!sample)
		{
			sample = std::make_shared<samples::sample_type>();
			sample->id = columns[header[sample_name_field]];

			std::vector<std::string> tags;
			boost::split(tags, columns[header[sample_tags_field]], boost::is_any_of(","));
			for(auto const& tag: tags)
			{
				std::vector<std::string> key_value;
				boost::split(key_value, tag, boost::is_any_of(":"));
				if(key_value.size() == 1)
				{
					sample->tags.insert(std::make_pair(key_value[0], ""));
				}
				else if(key_value.size() == 2)
				{
					sample->tags.insert(std::make_pair(key_value[0], key_value[1]));
				}
			}
		}

		auto& rearrangement(rearrangements[columns[header[rearrangement_field]]]);
		if(!rearrangement)
		{
			rearrangement = std::make_shared<rearrangements::rearrangement_type>();
			rearrangement->id = rearrangements.size();
			rearrangement->nucleotides = columns[header[rearrangement_field]];
		}

		auto sequence(std::make_shared<sequences::sequence_type>());
		sequence->id = ++ sequence_id;
		sequence->rearrangement = rearrangement;
		sequence->reads = boost::lexical_cast<std::size_t>(columns[header[reads_field]]);
		sequence->v_family = columns[header[v_family_field]];
		sequence->d_family = columns[header[d_family_field]];
		sequence->j_family = columns[header[j_family_field]];
		sample->sequences.push_back(sequence);
	}

	std::cout << "Processed " << rearrangements.size() << " unique rearrangements in " << sequence_id << " sequences" << std::endl;
}


}	// namespace immunoseq