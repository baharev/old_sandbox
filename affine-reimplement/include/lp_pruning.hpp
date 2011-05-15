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

#ifndef LP_PRUNING_HPP_
#define LP_PRUNING_HPP_

#include <vector>

namespace asol {

class lp_impl;

class lp_pruning {

public:

	lp_pruning(lp_impl* lp, const std::vector<int>& index_set);

private:

	lp_pruning(const lp_pruning& );
	lp_pruning& operator=(const lp_pruning& );

	enum decision { NO_CANDIDATE, LOWER_BND, UPPER_BND };

	void init_reverse_index_set();
	void init_bounds();
	void init_reduced_costs();

	void examine_col(int i);
	void examine_lb(int i, double val);
	void examine_ub(int i, double val);
	void mark_min_solved(int i);
	void mark_max_solved(int i);

	void dbg_examine_results() const;

	double col_val_with_check(int i);

	decision select_candidate();

	void count_solved() const;

	void solve_for_lb();

	void solve_for_ub();

	void prune();

	void prune_all();

	lp_impl* const lp;

	const std::vector<int>& index_set;

	const size_t size;

	std::vector<int> reverse_index_set;

	std::vector<char> min_solved;
	std::vector<char> max_solved;

	std::vector<double> lo;
	std::vector<double> up;

	std::vector<std::vector<double> > d_min;
	std::vector<std::vector<double> > d_max;

	int skipped;
	int index_min;
	int index_max;
	double closest_min;
	double closest_max;
};

}

#endif // LP_PRUNING_HPP_
