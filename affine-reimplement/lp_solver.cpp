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
#include <iostream>
#include <limits>
#include "lp_solver.hpp"
#include "affine.hpp"
#include "diagnostics.hpp"
#include "lp_impl.hpp"

using std::fabs;

namespace asol {

lp_solver::lp_solver() : lp(new lp_impl), N_VARS(-1), TINY(1.0e-7) {

}

void lp_solver::reset() {

	lp->reset();

	lp->add_cols(N_VARS);
}

void lp_solver::set_number_of_vars(int n) {

	ASSERT2( n>0 && N_VARS==-1, "n, N_VARS: " << n << ", " << N_VARS);

	// TODO Add as many columns?
	N_VARS = n;
}

void lp_solver::reset_col_arrays(int size) {

	col_index.clear();
	col_coeff.clear();

	ASSERT(size > 0);

	col_index.reserve(size);
	col_coeff.reserve(size);

	col_index.push_back(0);
	col_coeff.push_back(0.0);
}

int lp_solver::col_size() const {

	ASSERT(col_index.size()==col_coeff.size());

	return static_cast<int>(col_index.size());
}

void lp_solver::add_equality_constraint(const affine& x, const double value) {

	//std::cout << "x:\n" << x << std::endl;

	const row_info row = compute_row_info(x, value);

	const int n = x.size();

	reset_col_arrays(n);

	int i=1;

	for ( ; i<n; ++i) {

		const epsilon& e = x.noise_vars.at(i);

		if (e.index > N_VARS) {

			break;
		}

		if (std::fabs(e.coeff) > row.tiny) {

			col_index.push_back(e.index);
			col_coeff.push_back(e.coeff);
		}
	}

	ASSERT2(i >= col_size(), "i, size: "<<i<<", "<<col_size());

	lp->add_eq_row(&col_index.at(0), &col_coeff.at(0), col_size()-1, row.lb, row.ub);

	set_col_bounds();

	//lp->dump("lp_dump.txt");

	// TODO Call dual simplex here!
	lp->check_feasibility();
}

const lp_solver::row_info lp_solver::compute_row_info(const affine& x, const double value) const {

	const row_rad_max_aij row = get_row_rad_max_aij(x);

	const double mid = value - x.central_value();

	const double row_lb = mid - row.rad;

	const double row_ub = mid + row.rad;

	const double row_max = std::max(fabs(row_lb), fabs(row_ub));

	double max_aij = row.max_aij;

	if (row_max > max_aij) {

		max_aij = row_max;
	}

	ASSERT(max_aij > 0.0);

	return row_info(row_lb, row_ub, max_aij*TINY);
}

const lp_solver::row_rad_max_aij lp_solver::get_row_rad_max_aij(const affine& x) const {

	double max_aij = 0.0;

	double rad = 0.0;

	for (int i=1; i<x.size(); ++i) {

		const epsilon& e = x.noise_vars.at(i);

		const double aij = fabs(e.coeff);

		if (e.index > N_VARS) { // condense non-vars, needed for row bounds

			rad += aij;
		}
		else if (aij > max_aij) { // e.index <= N_VARS, find max a_ij

			max_aij = aij;
		}
	}

	return row_rad_max_aij(rad, max_aij);
}

// TODO Rename it to refresh?
void lp_solver::check_feasibility() {

	lp->check_feasibility();
}

void lp_solver::set_col_bounds() {

	const int n = col_size() - 1;

	for (int i=1; i<=n; ++i) {

		const int index = col_index.at(i);

		lp->set_col_bounds(index, -1, 1);
	}
}

void lp_solver::show_iteration_count() const {

	lp->show_iteration_count();
}

lp_solver::~lp_solver() {

	delete lp;
}

void lp_solver::free_environment() {

	lp_impl::free_environment();
}

}
