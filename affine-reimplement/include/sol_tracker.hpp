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

#ifndef SOL_TRACKER_HPP_
#define SOL_TRACKER_HPP_

#include <vector>
#include <utility>
#include "problem.hpp"

namespace asol {

class builder;
class interval;

class sol_tracker {

public:

	explicit sol_tracker(const problem<builder>* prob);

	explicit sol_tracker(const DoubleArray2D& solutions);

	void save_containment_info(const std::vector<interval>* v);

	void print_containment_statistics() const;

	bool contains_solution() const;

	void dump_previous_v() const;

	void check_transitions_since_last_call(const std::vector<interval>* v);

	~sol_tracker();

private:

	enum type { NOT, EASY, STRICT};

	typedef std::pair<type,int> containment_info;

	typedef std::vector<std::vector<double> >::const_iterator const_itr;

	void save_containment(const const_itr begin, const const_itr end);

	const containment_info check_sol() const;

	int first_not_strictly_contained() const;

	int first_not_easily_contained(const int from) const;

	void check_all_sol(const const_itr begin, const const_itr end);

	void check_transition(const containment_info status, int sol_index) const;

	void show_component(const char* msg, const int index) const;

	const int n_vars;

	DoubleArray2D solutions;

	containment_info STRICT_CONTAINMENT;

	std::vector<containment_info> containment;

	const std::vector<interval>* v;

	std::vector<double>::const_iterator sol;

	std::vector<interval> previous_v;

};

}

#endif // SOL_TRACKER_HPP_
