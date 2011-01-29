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

//#include <iostream>
//#include <iterator>
#include "combination.hpp"
#include "diagnostics.hpp"

namespace asol {

combination::combination() {

	counter_max = size = position = -1;
}

combination::combination(int index_size, int parts_to_generate) {

	ASSERT(index_size>0);

	size = index_size;

	counter_max = parts_to_generate-1;

	counter.resize(size, 0);

	counter.at(0) = -1;

	position = 0;
}

bool combination::step_counters() {

	ASSERT(size>0);

	bool overflow = false;

	while (has_more_counters() && counter_at_max()) {

		++position;

		overflow = true;
	}

	if (!has_more_counters()) {

		return false;
	}

	next(overflow);

	return true;
}

bool combination::has_more_counters() const {

	return position < size;
}

bool combination::counter_at_max() const {

	return counter.at(position) == counter_max;
}

void combination::next(const bool overflow) {

	if (overflow) {

		handle_overflow();
	}
	else {

		++counter.at(position);
	}
/*
	using namespace std;

	copy(counter.begin(), counter.end(), ostream_iterator<int>(cout, "\t"));

	cout << endl;
*/
	return;
}

void combination::handle_overflow() {

	for (int i=0; i<position; ++i) {

		counter.at(i)=0;
	}

	++counter.at(position);

	position = 0;
}

const std::vector<int>& combination::counters() const {

	return counter;
}

}
