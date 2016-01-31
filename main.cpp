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
		
		for(auto const& patient: patients)
		{
			std::cout << "Patient " << patient.second->id << ":" << std::endl;
			for(auto const& sample: patient.second->samples)
			{
				std::cout << "\t" << sample.second->sequences.size() << " sequences from sample " << sample.second->id << " dated " << sample.second->timestamp.date() << std::endl;
			}
		}

		return 0;
	}

	catch(...)
	{
		std::cout << boost::current_exception_diagnostic_information() << std::endl;
		return 1;
	}
}