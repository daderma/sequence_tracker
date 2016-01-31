#include "edit_distance.hpp"
#include "samples.hpp"
#include "windows.hpp"
#include "immunoseq.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include <iostream>
#include <set>


typedef std::multimap<std::pair<std::size_t, std::int64_t>, sequences::sequence_ptr_type> distances_type;
std::size_t const maximum_distance(2);


void next_iteration(distances_type const& previous, sequences::sequences_type const& current, distances_type& next)
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
			next.insert(std::make_pair(std::make_pair(i.second, i.first->id), i.first));
		}
	}
}


typedef std::vector<samples::sample_ptr_type> stack_type;
void process(samples::samples_type const& samples, samples::sample_ptr_type const& parent, distances_type const& previous, std::string const& report, boost::filesystem::path const& directory, stack_type& stack)
{
	for(auto const& current: samples)
	{
		std::set<samples::sample_ptr_type> parents;
		for(auto const& tag: current.second->tags)
		{
			if(boost::iequals(tag.first, "ParentSample"))
			{
				auto const iter(samples.find(tag.second));
				if(iter == samples.end())
				{
					throw std::runtime_error("Sample " + current.first + " references non-existant parent " + tag.second);
				}
				else
				{
					parents.insert(iter->second);
				}
			}
		}

		if(!parent && !parents.empty())
		{
			continue;
		}

		if(parent && !parents.count(parent))
		{
			continue;
		}

		stack.push_back(current.second);
		for(std::size_t i(0); i < stack.size(); ++ i)
		{
			std::cout << "\t";
		}
		std::cout << current.first << ": " << current.second->sequences.size() << " sequences";
		distances_type next;
		next_iteration(previous, current.second->sequences, next);
		std::cout << " -> " << next.size() << " in next iteration" << std::endl;
		
		std::ostringstream buffer;
		buffer << report;
		buffer << current.second->sequences.size() << " sequences from sample " << current.first << std::endl;
		for(auto const& distance: next)
		{
			buffer << "\t" << distance.first.first << "\t" << distance.second << std::endl;
		}

		std::string filename;
		for(auto const& iteration: stack)
		{
			if(filename.empty())
			{
				filename = iteration->id;
			}
			else
			{
				filename += "_" + iteration->id;
			}
		}
		boost::filesystem::ofstream stream(directory / (filename + ".txt"), std::ios::trunc);
		stream << buffer.str();
		process(samples, current.second, next, buffer.str(), directory, stack);
		stack.pop_back();
	}
}


void try_directory(boost::filesystem::path const& directory)
{
	if(boost::filesystem::is_directory(directory))
	{
		samples::samples_type samples;

		boost::filesystem::directory_iterator end;
		for(boost::filesystem::directory_iterator iter(directory); iter != end; ++ iter)
		{
			if(boost::filesystem::is_regular_file(iter->path()) && boost::iequals(iter->path().extension().string(), ".tsv"))
			{
				immunoseq::load(iter->path(), samples);	
			}
		}

		std::cout << "Performing iterative analysis:" << std::endl;
		distances_type previous;
		std::string report;
		stack_type stack;
		process(samples, samples::sample_ptr_type(), previous, report, directory, stack);
	}
	else
	{
		std::cout << "Not a directory: " << directory << std::endl;
	}
}


int main(int argc, char* argv[])
{
	try
	{
		if(argc > 1)
		{
			for(int i(1); i < argc; ++ i)
			{
				try_directory(argv[i]);
			}
		}
		else
		{
			std::cout << "This program processes a collection of immunoSEQ samples located in tab-separated-values (.tsv) files. Example:" << std::endl;
			std::cout << "samples" << std::endl << "\tproject1.tsv" << std::endl << "\tproject2.tsv" << std::endl << "\tproject3.tsv" << std::endl << "\t..." << std::endl << std::endl;
			std::cout << "You can pass the directory (i.e. \"samples\" above) containing the files you wish to process on the command line." << std::endl;
			std::cout << "For now you can manually select a directory to process..." << std::endl;
			try_directory(windows::select_directory());
		}
		
		std::cout << std::endl << "Press (almost) any key to continue." << std::endl;
		std::cin.get();
		return 0;
	}

	catch(...)
	{
		std::cout << boost::current_exception_diagnostic_information() << std::endl;
		return 1;
	}
}