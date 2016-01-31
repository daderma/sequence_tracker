#include "edit_distance.hpp"
#include "patients.hpp"
#include "windows.hpp"
#include "immunoseq.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include <iostream>


void try_directory(boost::filesystem::path const& directory, patients::patients_type& patients)
{
	if(boost::filesystem::is_directory(directory))
	{
		boost::filesystem::directory_iterator end;
		for(boost::filesystem::directory_iterator iter(directory); iter != end; ++ iter)
		{
			if(boost::filesystem::is_regular_file(iter->path()) && boost::iequals(iter->path().extension().string(), ".tsv"))
			{
				immunoseq::load(iter->path(), patients);	
			}
		}
	}
	else
	{
		std::cout << "Not a directory: " << directory << std::endl;
	}
}


void top_clone_distances(std::ostream& stream, sequences::sequences_type const& sequences, std::size_t const& limit)
{
	sequences::sequence_ptr_type top_clone;
	for(auto const& sequence: sequences)
	{
		if(!top_clone || top_clone->reads < sequence->reads)
		{
			top_clone = sequence;
		}
	}
	
	std::multimap<std::size_t, sequences::sequence_ptr_type> distances;
	for(auto const& sequence: sequences)
	{
		if(sequence != top_clone)
		{
			distances.insert(std::make_pair(edit_distance::levenshtein(top_clone->rearrangement, sequence->rearrangement), sequence));
		}
	}

	stream << "\t\tTop\t" << top_clone->rearrangement << " (" << top_clone->reads << " reads)" << std::endl;

	std::size_t result(0);
	for(auto const& distance: distances)
	{
		if(++ result < limit)
		{
			stream << "\t\t" << distance.first << "\t" << distance.second->rearrangement << " (" << distance.second->reads << " reads)" << std::endl;
		}
		else
		{
			break;
		}
	}
}


int main(int argc, char* argv[])
{
	try
	{
		patients::patients_type patients;
		if(argc > 1)
		{
			for(int i(1); i < argc; ++ i)
			{
				try_directory(argv[i], patients);
			}
		}
		else
		{
			std::cout << "This program processes a collection of immunoSEQ samples located in tab-separated-values (.tsv) files. Example:" << std::endl;
			std::cout << "samples" << std::endl << "\tproject1.tsv" << std::endl << "\tproject2.tsv" << std::endl << "\tproject3.tsv" << std::endl << "\t..." << std::endl << std::endl;
			std::cout << "You can pass the directory (i.e. \"samples\" above) containing the files you wish to process on the command line." << std::endl;
			std::cout << "For now you can manually select a directory to process..." << std::endl;
			try_directory(windows::select_directory(), patients);
		}
		

		std::cout << "Creating report";

		std::ofstream output("top-clone.txt", std::ios::trunc);
		for(auto const& patient: patients)
		{
			output << "Patient " << patient.second->id << ":" << std::endl;
			
			sequences::sequences_type aggregated;
			for(auto const& sample: patient.second->samples)
			{
				std::cout << ".";

				output << "\t" << sample.second->sequences.size() << " sequences from sample " << sample.second->id << " dated " << sample.second->timestamp.date() << std::endl;
				top_clone_distances(output, sample.second->sequences, 20);
				aggregated.insert(aggregated.end(), sample.second->sequences.begin(), sample.second->sequences.end());
			}

			output << "\t" << aggregated.size() << " sequences from all patient samples" << std::endl;
			top_clone_distances(output, aggregated, 20);
		}

		std::cout << std::endl;

		return 0;
	}

	catch(...)
	{
		std::cout << boost::current_exception_diagnostic_information() << std::endl;
		return 1;
	}
}