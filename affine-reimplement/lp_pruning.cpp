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

#include <iostream>
#include <limits>
#include "lp_pruning.hpp"
#include "lp_impl.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"

using std::vector;

namespace {

const double TOL_LP_PRUNING_SOLVED = 1.0e-8;

}

namespace asol {

lp_pruning::lp_pruning(lp_impl* lp , const vector<int>& index_set)
:
lp(lp),
index_set(index_set),

size(index_set.size()),

min_solved(size, 'n'),
max_solved(size, 'n'),

lo(size, 10.0),
up(size,-10.0),

skipped(0)

{
	init_reverse_index_set();

	init_bounds();

	init_reduced_costs();

	prune_all();
}

void lp_pruning::init_reverse_index_set() {

	const int n_cols = lp->num_cols();

	reverse_index_set.resize(n_cols+1, -1);

	for (size_t i=0; i<size; ++i) {

		const int index = index_set.at(i);

		reverse_index_set.at(index) = i;
	}
}

void lp_pruning::init_bounds() {

	for (size_t i=0; i<size; ++i) {

		const int index = index_set.at(i);

		const double lb = lp->col_lb(index);
		const double ub = lp->col_ub(index);

		ASSERT2(lb < ub, "lb, ub: "<<lb<<", "<<ub);

		lo.at(i) = lb;
		up.at(i) = ub;
	}
}

void lp_pruning::init_reduced_costs() {

	const int n_cols = lp->num_cols();
	const int n_rows = lp->num_rows();

	d_min.resize(1+n_rows);
	d_max.resize(1+n_rows);

	for (int i=0; i<=n_rows; ++i) {

		d_min.at(i) = vector<double>(1+n_cols, 0.0);
		d_max.at(i) = vector<double>(1+n_cols, 0.0);
	}
}

void lp_pruning::prune_all() {

	// TODO mark_narrow_solved();

	try {

		prune();
	}
	catch (infeasible_problem& ) {

		throw numerical_problems();
	}

	count_solved();
}

void lp_pruning::prune() {

	count_solved();

	size_t lp_call = 0;

	subroblem next;

	while ( (next=select_candidate()) != NO_MORE ) {

		count_solved();

		if (next == MIN_SUBPROBLEM) {

			solve_for_lb();
		}
		else {

			solve_for_ub();
		}

		++lp_call;
	}

	ASSERT(lp_call+skipped==2*size);
}

void lp_pruning::count_solved() const {

	int solved = 0;

	for (size_t i=0; i<size; ++i) {

		if (min_solved.at(i)=='y') {

			++solved;
		}

		if (max_solved.at(i)=='y') {

			++solved;
		}
	}

	std::cout << "Solved: " << solved << "/" << 2*size << ", ";
	std::cout << "skipped: " << skipped << std::endl;
}

lp_pruning::subroblem lp_pruning::select_candidate() {

	closest_min = closest_max = std::numeric_limits<double>::max();

	index_min = index_max = -1;

	for (size_t i=0; i<size; ++i) {

		examine_lb(i);

		examine_ub(i);
	}

	dbg_selection_results();

	subroblem next;

	if (closest_min < closest_max) {

		next = MIN_SUBPROBLEM;
	}
	else if (index_max!=-1) { // closest_max <= closest_min

		next = MAX_SUBPROBLEM;
	}
	else {

		next = NO_MORE;
	}

	return next;
}

void lp_pruning::examine_lb(int i) {

	if (min_solved.at(i)=='y') {

		return;
	}

	double val = lp->col_val(index_set.at(i));

	double distance_from_lb = val-lo.at(i);

	ASSERT(distance_from_lb >= 0);

	if (distance_from_lb <= TOL_LP_PRUNING_SOLVED) {

		min_solved.at(i) = 'y';

		++skipped;
	}
	else if (distance_from_lb < closest_min) {

		closest_min = distance_from_lb;

		index_min = i;
	}
}

void lp_pruning::examine_ub(int i) {

	if (max_solved.at(i)=='y') {

		return;
	}

	double val = lp->col_val(index_set.at(i));

	double distance_from_ub = up.at(i)-val;

	ASSERT(distance_from_ub >= 0);

	if (distance_from_ub <= TOL_LP_PRUNING_SOLVED) {

		max_solved.at(i) = 'y';

		++skipped;
	}
	else if (distance_from_ub < closest_max) {

		closest_max = distance_from_ub;

		index_max = i;
	}
}

void lp_pruning::solve_for_lb() {

	min_solved.at(index_min) = 'y';

	lp->tighten_col_lb(index_set.at(index_min), lo.at(index_min));
}

void lp_pruning::solve_for_ub() {

	max_solved.at(index_max) = 'y';

	lp->tighten_col_ub(index_set.at(index_max), up.at(index_max));
}

void lp_pruning::dbg_selection_results() const {

	using namespace std;

	cout << "min: ";
	if (index_min!=-1) {
		cout << index_set.at(index_min) << ", val: " << closest_min << endl;
	}
	else {
		ASSERT(closest_min==numeric_limits<double>::max());
		cout << " (no more)" << endl;
	}

	cout << "max: ";
	if (index_max!=-1) {
		cout << index_set.at(index_max) << ", val: " << closest_max << endl;
	}
	else {
		ASSERT(closest_max==numeric_limits<double>::max());
		cout << " (no more)" << endl;
	}
}

}
