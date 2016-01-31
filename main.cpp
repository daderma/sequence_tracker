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


typedef std::multimap<std::size_t, sequences::sequence_ptr_type> distances_type;


void next_iteration(distances_type const& previous, sequences::sequences_type const& current, std::size_t const& maximum_distance, distances_type& next)
{
	std::map<sequences::sequence_ptr_type, std::size_t> inverse;

	if(previous.empty())
	{
		sequences::sequence_ptr_type top;
		for(auto const& c: current)
		{
			if(!top || top->reads < c->reads)
			{
				top = c;
			}
		}

		for(auto const& c: current)
		{
			inverse[c] = edit_distance::levenshtein(top->rearrangement->nucleotides, c->rearrangement->nucleotides);
		}
	}
	else
	{
		for(auto const& p: previous)
		{
			for(auto const& c: current)
			{
				if(inverse.count(c))
				{
					inverse[c] = std::min(inverse[c], edit_distance::levenshtein( p.second->rearrangement->nucleotides, c->rearrangement->nucleotides));
				}
				else
				{
					inverse[c] = edit_distance::levenshtein( p.second->rearrangement->nucleotides, c->rearrangement->nucleotides);
				}
			}
		}
	}

	for(auto const& i: inverse)
	{
		if(i.second <= maximum_distance)
		{
			next.insert(std::make_pair(i.second, i.first));
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

		std::ofstream report("iterations.txt", std::ios::trunc);
		for(auto const& patient: patients)
		{
			report << "Patient " << patient.second->id << ":" << std::endl;

			distances_type previous;
			for(auto const& sample: patient.second->samples)
			{
				std::cout << ".";
				
				report << "\t" << sample.second->sequences.size() << " sequences from sample " << sample.second->id << " dated " << sample.second->timestamp.date() << std::endl;

				distances_type next;
				next_iteration(previous, sample.second->sequences, 3, next);
				
				for(auto const& distance: next)
				{
					report << "\t\t" << distance.first << "\t" << distance.second << std::endl;
				}
	
				previous = next;
			}
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