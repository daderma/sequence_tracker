#include "edit_distance.hpp"
#include <algorithm>
#include <numeric>


namespace edit_distance
{


std::size_t levenshtein(std::string const& left, std::string const& right)
{
	std::size_t const left_size(left.size());
	std::size_t const right_size(right.size());
	std::size_t column_start(1);
	auto column(new std::size_t[left_size + 1]);
	std::iota(column + column_start, column + left_size + 1, column_start);
	
	for(auto x(column_start); x <= right_size; ++ x)
	{
		column[0] = x;
		auto last_diagonal(x - column_start);
		for(auto y(column_start); y <= left_size; ++ y)
		{
			auto old_diagonal(column[y]);
			auto possibilities = 
			{
				column[y] + 1,
				column[y - 1] + 1,
				last_diagonal + (left[y - 1] == right[x - 1] ? 0 : 1)
			};
			column[y] = std::min(possibilities);
			last_diagonal = old_diagonal;
		}
	}
	auto result = column[left_size];
	delete[] column;
	return result;
}


}	// namespace edit_distance