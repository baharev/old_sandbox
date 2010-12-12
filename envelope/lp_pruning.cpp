//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010 Ali Baharev
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
#include <algorithm>
#include <assert.h>
#include "lp_pruning.hpp"
#include "lp_impl.hpp"
#include "exceptions.hpp"
#include "constants.hpp"

namespace asol {

lp_pruning::lp_pruning(lp_impl* lpmin, lp_impl* lpmax, int to_index)
: n(to_index),
  min_solved(new bool[1+n]),
  max_solved(new bool[1+n]),
  lo(new double[1+n]),
  up(new double[1+n])
{
	assert(n<=lpmax->n_cols());
	lp_min = lpmin;
	lp_max = lpmax;
	skipped = 0;
}

lp_pruning::~lp_pruning() {

	delete[] min_solved;
	delete[] max_solved;
	delete[] lo;
	delete[] up;
}

void lp_pruning::init() {

	lp_min->refresh_basis();
	lp_max->refresh_basis();

	for (int i=1; i<=n; ++i) {

		min_solved[i] = max_solved[i] = false;

		const double lb = lp_min->col_lb(i);
		const double ub = lp_min->col_ub(i);

		assert(lb==lp_max->col_lb(i));
		assert(ub==lp_max->col_ub(i));
		assert(lb <= ub);

		lo[i] = lb;
		up[i] = ub;
	}
}

void lp_pruning::mark_narrow_solved() {

	for (int i=1; i<=n; ++i) {

		if (too_narrow(lo[i], up[i])) {

			min_solved[i] = max_solved[i] = true;
		}
	}
}

void lp_pruning::examine_col(int i) {

	double diam = up[i]-lo[i];

	if (diam==0 || (min_solved[i] && max_solved[i])) {

		return;
	}

	double val_min = lp_min->get_col_val(i);

	double val_max = lp_max->get_col_val(i);

	double threshold = TOL_PRUNING_SOLVED*diam;

	double offcenter_lb = std::min(val_min-lo[i], val_max-lo[i])/diam;

	examine_lb(i, offcenter_lb, threshold);

	double offcenter_ub = std::min(up[i]-val_min, up[i]-val_max)/diam;

	examine_ub(i, offcenter_ub, threshold);
}

void lp_pruning::examine_lb(int i, double offcenter_lb, double threshold) {

	if (!min_solved[i] && (offcenter_lb < threshold)) {

		min_solved[i] = true;

		++skipped;
	}

	if (!min_solved[i] && offcenter_lb < closest_min) {

		closest_min = offcenter_lb;

		index_min = i;
	}
}

void lp_pruning::examine_ub(int i, double offcenter_ub, double threshold) {

	if (!max_solved[i] && (offcenter_ub < threshold)) {

		max_solved[i] = true;

		++skipped;
	}

	if (!max_solved[i] && offcenter_ub < closest_max) {

		closest_max = offcenter_ub;

		index_max = i;
	}
}

int lp_pruning::select_candidate() {

	closest_min = closest_max = std::numeric_limits<double>::max();

	index_min = index_max = -1;

	for (int i=1; i<=n; ++i) {

		examine_col(i);
	}

	//std::cout << "min: " << index_min << ", val: " << closest_min << std::endl;
	//std::cout << "max: " << index_max << ", val: " << closest_max << std::endl;

	decision ret_val;

	if (closest_min < closest_max) {

		ret_val = LOWER_BND;
	}
	else if (index_max!=-1){

		ret_val = UPPER_BND;
	}
	else {

		ret_val = NO_CANDIDATE;
	}

	return ret_val;
}

void lp_pruning::count_solved() const {

	int solved = 0;

	for (int i=1; i<=n; ++i) {

		if (min_solved[i]) {

			++solved;
		}

		if (max_solved[i]) {

			++solved;
		}
	}

	//std::cout << "Solved: " << solved << "/" << 2*n << ", ";
	//std::cout << "skipped: " << skipped << std::endl;
}

void lp_pruning::solve_for_lb() {

	min_solved[index_min] = true;

	double value_in_min = lp_min->get_col_val(index_min);
	double value_in_max = lp_max->get_col_val(index_min);

	lp_impl* const lp = (value_in_min < value_in_max)? lp_min : lp_max;

	lp->tighten_col_lb(index_min, lo[index_min]);
}

void lp_pruning::solve_for_ub() {

	max_solved[index_max] = true;

	double value_in_min = lp_min->get_col_val(index_max);
	double value_in_max = lp_max->get_col_val(index_max);

	lp_impl* const lp = (value_in_min < value_in_max)? lp_max : lp_min;

	lp->tighten_col_ub(index_max, up[index_max]);
}

void lp_pruning::prune() {

	count_solved();

	int candidate;

	while ( (candidate=select_candidate()) != NO_CANDIDATE ) {

		count_solved();

		if (candidate == LOWER_BND) {

			solve_for_lb();
		}
		else {

			solve_for_ub();
		}
	}
}

void lp_pruning::prune_all() {

	init();

	mark_narrow_solved();

	try {

		prune();
	}
	catch (asol::infeasible_problem& ) {
		std::cout << "Warning: numerical problems " << __FILE__ << " ";
		std::cout << __LINE__ << std::endl;
		throw asol::numerical_problems();
	}

	count_solved();

	// TODO Check consistency
	// TODO Check if progress is sufficient

}

}
