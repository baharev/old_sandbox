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

#ifndef LP_PRUNING_HPP_
#define LP_PRUNING_HPP_

namespace asol {

	class interval;
}

namespace lp_solver {

class lp_impl;

class lp_pruning {

public:

	lp_pruning(lp_impl* lpmin, lp_impl* lpmax, int to_index);

	void prune_all();

	friend void copy_bounds(const lp_pruning& lp, asol::interval* bounds);

	~lp_pruning();

private:

	lp_pruning(const lp_pruning& );
	lp_pruning& operator=(const lp_pruning& );

	void init();
	void mark_narrow_solved();
	void count_solved() const;
	void examine_col(int i);
	void examine_lb(int i, double offcenter_lb, double threshold);
	void examine_ub(int i, double offcenter_ub, double threshold);
	int  select_candidate();
	void solve_for_lb();
	void solve_for_ub();

	const int n;
	bool* const min_solved;
	bool* const max_solved;
	double* const lo;
	double* const up;
	lp_impl* lp_min;
	lp_impl* lp_max;

	double closest_min;
	double closest_max;
	int index_min;
	int index_max;

	enum decision { NO_CANDIDATE, LOWER_BND, UPPER_BND };
};

}

#endif /* LP_PRUNING_HPP_ */
