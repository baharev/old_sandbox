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

#include <cmath>
#include <limits>
#include "lp_solver.hpp"
#include "affine.hpp"
#include "diagnostics.hpp"
#include "lp_impl.hpp"

using std::fabs;

namespace asol {

lp_solver::lp_solver() : lp(new lp_impl), N_VARS(-1), TINY(1.0e-7) {

}

void lp_solver::set_number_of_vars(int n) {

	ASSERT2( n>0 && N_VARS==-1, "n, N_VARS: " << n << ", " << N_VARS);

	// TODO Add as many columns?
	N_VARS = n;
}

void lp_solver::reset_col_arrays() {

	col_index.clear();
	col_coeff.clear();

	col_index.push_back(0);
	col_coeff.push_back(0.0);
}

void lp_solver::reserve_col_arrays(int size) {

	ASSERT(size > 0);

	col_index.reserve(size);
	col_coeff.reserve(size);
}

void lp_solver::add_equality_constraint(const affine& x, const double value) {

	const int n = x.size();

	double max_aij = 0.0;

	int i=1;

	for ( ; i<n; ++i) {

		const epsilon& e = x.noise_vars.at(i);

		if (e.index > N_VARS) {

			break;
		}

		const double aij = fabs(e.coeff);

		if (aij > max_aij) {

			max_aij = aij;
		}
	}

	double rad = 0.0;

	for ( ; i<n; ++i) {

		const epsilon& e = x.noise_vars.at(i);

		const double aij = fabs(e.coeff);

		rad += aij;
	}

	const double mid = x.central_value();

	const double row_lb = -mid - rad;

	const double row_ub = -mid + rad;

	const double row_max = std::max(fabs(row_lb), fabs(row_ub));

	if (row_max > max_aij) {

		max_aij = row_max;
	}

	ASSERT(max_aij > 0);

	reset_col_arrays();

	reserve_col_arrays(n);
}

lp_solver::~lp_solver() {

	delete lp;
}

void lp_solver::free_environment() {

	lp_impl::free_environment();
}

}
