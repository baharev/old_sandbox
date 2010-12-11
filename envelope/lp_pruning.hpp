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

namespace lp_solver {

class lp_impl;

class lp_pruning {

public:

	lp_pruning(lp_impl* lp_min, lp_impl* lp_max);

	void prune_all();

	~lp_pruning();

private:

	lp_pruning(const lp_pruning& );
	lp_pruning& operator=(const lp_pruning& );

	void init();
	void mark_narrow_solved();

	const int n;
	bool* const min_solved;
	bool* const max_solved;
	double* const lo;
	double* const up;
	lp_impl* lp_min;
	lp_impl* lp_max;
};

}

#endif /* LP_PRUNING_HPP_ */
