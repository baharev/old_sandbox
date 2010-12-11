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

#include <assert.h>
#include "lp_pruning.hpp"
#include "lp_impl.hpp"

namespace lp_solver {

lp_pruning::lp_pruning(lp_impl* lp_min, lp_impl* lp_max)
: n(lp_min->n_cols()),
  min_solved(new bool[1+n]),
  max_solved(new bool[1+n]),
  lo(new double[1+n]),
  up(new double[1+n])
{
	assert(n==lp_max->n_cols());
	this->lp_min = lp_min;
	this->lp_max = lp_max;
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

void lp_pruning::prune_all() {

	init();

	mark_narrow_solved();
}

}
