//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
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
#include <iostream>
#include "sol_tracker.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "interval.hpp"
#include "vector_dump.hpp"

using namespace std;

namespace asol {

sol_tracker::sol_tracker(const problem<builder>* prob)
: n_vars(prob->number_of_variables()),
  solutions(prob->solutions())
{
	// TODO Make number_of_stored_solutions() return size()
	const int size = static_cast<int>(solutions.size());
	ASSERT2(size==prob->number_of_stored_solutions(),"size: "<<size);

	STRICT_CONTAINMENT = make_pair(STRICT, -1); // FIXME Why fails in the init list?
	v = 0;
}

sol_tracker::sol_tracker(const DoubleArray2D& sols)
: n_vars(sols.at(0).size()), solutions(sols)
{
	STRICT_CONTAINMENT = make_pair(STRICT, -1);
	v = 0;
}

sol_tracker::~sol_tracker() {
	// Just to make the compiler shut-up
}

void sol_tracker::save_containment_info(const std::vector<interval>* current_v) {

	ASSERT2(!solutions.empty(),"set solutions first");

	v = current_v;

	save_containment(solutions.begin(), solutions.end());

	previous_v.assign(v->begin(), v->end());

	v = 0;

	print_containment_statistics();
}

void sol_tracker::print_containment_statistics() const {

	int strict = count(containment.begin(), containment.end(), STRICT_CONTAINMENT);

	cout << "Strictly contains " << strict << " of " << containment.size();
	cout << " solutions";

	if ( strict == 1) {
		int pos = find(containment.begin(), containment.end(), STRICT_CONTAINMENT)
				     - containment.begin();
		cout << " (" << pos+1 << ")";
	}

	cout << endl;
}

bool sol_tracker::contains_solution() const {

	return find(containment.begin(),containment.end(),STRICT_CONTAINMENT)!=containment.end();
}

void sol_tracker::save_containment(const const_itr begin, const const_itr end) {

	ASSERT(containment.empty());

	for (const_itr i = begin; i != end; ++i) {

		sol = i->begin();

		containment_info status = check_sol();

		containment.push_back( status );
	}
}

const sol_tracker::containment_info sol_tracker::check_sol() const {

	const int i = first_not_strictly_contained();

	if (i==n_vars) {

		return STRICT_CONTAINMENT;
	}

	const int j = first_not_easily_contained(i);

	if (j==n_vars) {

		return containment_info(EASY, i);
	}
	else {

		return containment_info(NOT, j);
	}
}

int sol_tracker::first_not_strictly_contained() const {

	for (int i=0; i<n_vars; ++i) {

		if (! v->at(i).contains(sol[i]) ) {

			return i;
		}
	}

	return n_vars;
}

int sol_tracker::first_not_easily_contained(const int from) const {

	for (int i=from ; i<n_vars; ++i) {

		if (!easy_containment(sol[i], v->at(i)) ) {

			return i;
		}
	}

	return n_vars;
}

void sol_tracker::check_transitions_since_last_call(const std::vector<interval>* current_v) {

	ASSERT2(!containment.empty(),"save containment info first");

	v = current_v;

	check_all_sol(solutions.begin(), solutions.end());

	previous_v.assign(v->begin(), v->end());

	v = 0;
}

void sol_tracker::check_all_sol(const const_itr begin, const const_itr end) {

	for (const_itr i = begin; i != end; ++i) {

		sol = i->begin();

		containment_info status = check_sol();

		check_transition(status, i-begin);

		containment.at(i-begin) = status;
	}
}

void sol_tracker::check_transition(const containment_info status, int sol_index) const {

	const int previous = containment.at(sol_index).first;

	const int current  = status.first;

	const int index = status.second;

	if (previous == STRICT && current == EASY) {

		show_component("Warning: strict containment became easy", index);
	}
	else if (previous == STRICT && current == NOT) {

		show_component("Error: solution lost, generating hexadecimal dump", index);
		// TODO Hexa dump then exit
		ASSERT(false);
	}
	else if (previous < current) {

		show_component("Error: the box grew, generating hexadecimal dump", index);
		// TODO Hexa dump then exit
		ASSERT(false);
	}
	else if (previous == current) {

		; // Allowed
	}
	else {

		ASSERT(false);
	}
}

void sol_tracker::show_component(const char* msg, const int index) const {

	ios_base::fmtflags flag_old = cout.setf(ios_base::scientific);

	streamsize         prec_old = cout.precision(16);

	cout << msg << endl;

	cout << index << ": " << previous_v.at(index) << endl;

	cout << index << ": " <<           v->at(index) << "  ";

	cout << sol[index] << endl;

	cout.flags(flag_old);

	cout.precision(prec_old);
}

void sol_tracker::dump_previous_v() const {

	ASSERT(!previous_v.empty());

	asol::dump(previous_v);
}

}
