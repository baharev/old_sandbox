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
#include <iostream>
#include "box_generator.hpp"
#include "combination.hpp"
#include "diagnostics.hpp"

using namespace std;

namespace asol {

box_generator::box_generator(IVector& vec, const IntVector& index_set, int equal_parts) :
v(vec), parts_to_generate(equal_parts)
{
	ASSERT2(equal_parts>=2,"minimum 2 parts should be generated, asked for "<<equal_parts);

	reserve(index_set.size());

	for_each(index_set.begin(), index_set.end(), bind1st(mem_fun(&box_generator::generate_parts), this));

	if (!index.empty()) {

		index_generator.reset(new combination(index.size(), equal_parts));
	}

	ASSERT(index.size()==parts.size());
}

box_generator::~box_generator() {
	// Do NOT remove, needed to generate dtor of auto_ptr
}

void box_generator::reserve(int index_set_size) {

	ASSERT(index_set_size>0);

	int size = parts_to_generate;

	while (--index_set_size) {

		size*=parts_to_generate;
	}

	index.reserve(size);

	parts.reserve(size);
}

void box_generator::generate_parts(int i) {

	const interval range = v.at(i);

	if (!is_narrow(range)) {

		index.push_back(i);

		cut_into_equal_parts(range.inf(), range.sup());
	}
}

void box_generator::cut_into_equal_parts(const double LB, const double UB) {

	const int n = parts_to_generate;

	IVector tmp_parts(n);

	double lo = LB;

	for (int k=1; k<n; ++k) {

		const double up = ((n-k)*LB+k*UB)/n;

		tmp_parts.at(k-1) = interval(lo, up);

		lo = up;
	}

	tmp_parts.at(n-1) = interval(lo, UB); // k == n

	parts.push_back(tmp_parts);
}

bool box_generator::empty() const {

	return parts.empty();
}

bool box_generator::set_next() {

	ASSERT(!empty());

	bool has_more = index_generator->step_counters();

	if (!has_more) {

		return false;
	}

	// TODO Set box here

	const IntVector& counters = index_generator->counters();

	for (int i=0; i<counters.size(); ++i) {

		cout << parts.at(i).at(counters.at(i)) << '\t';
	}

	cout << endl;

	return true;
}

}

