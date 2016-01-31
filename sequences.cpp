#include "sequences.hpp"


namespace sequences
{


std::ostream& operator<<(std::ostream& stream, sequence_ptr_type const& sequence)
{
	stream 
		<< sequence->rearrangement->nucleotides
		<< "  id " << sequence->rearrangement->id
		<< "; reads " << sequence->reads << 
		"; families " << sequence->v_family << " " << sequence->d_family << " " << sequence->j_family;
	return stream;
}


}	// namespace sequences