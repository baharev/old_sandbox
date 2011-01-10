//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010, 2011 Ali Baharev
// All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//==============================================================================

#include <algorithm>
#include <functional>
#include "combination.hpp"
#include "interval.hpp"

using namespace std;

namespace asol {

combination::combination(const Vector& index_bound, int equal_parts) :
index(   new int[index_bound.size()]),
max_part(new int[index_bound.size()]),
counter( new int[index_bound.size()]),
part   ( new interval[index_bound.size()]),
length(0)
{
	for_each(index_bound.begin(), index_bound.end(), bind1st(mem_fun(&combination::copy_if_not_narrow),this));
}

combination::~combination() {

	delete[] index;
	delete[] max_part;
	delete[] counter;
	delete[] part;
}

// TODO Implement it for intervals
bool is_narrow(const interval& ) {

	return false;
}

void combination::copy_if_not_narrow(const index_range ir) {

	if (is_narrow(ir.range())) {

		return;
	}

	index[length] = ir.index();
	part [length] = ir.range();

	++length;
}

}
