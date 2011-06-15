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

#ifndef LP_SOLVER_HPP_
#define LP_SOLVER_HPP_

#include <vector>

namespace asol {

class affine;
class lp_impl;

class lp_solver {

public:

	lp_solver();

	~lp_solver();

	void reset();

	void add_equality_constraint(const affine& x, const double value);

	void prune_upcoming_variables();

	void check_feasibility();

	void set_number_of_vars(int n);

	void set_affine_vars(std::vector<affine>* v);

	void set_pruning_indices(const std::vector<std::vector<int> >& indices_to_prune_after_constraint);

	// returns zero based index to be split, negative if none selected
	int prune(const std::vector<int>& index_set);

	void show_iteration_count() const;

	static void free_environment();

private:

	lp_solver(const lp_solver& );
	lp_solver& operator=(const lp_solver& );

	void make_index_set_one_based();

	struct row_info {
		row_info(double lb, double ub, double tiny) : lb(lb), ub(ub), tiny(tiny) { }
		double lb;
		double ub;
		double tiny;
	};

	struct row_rad_max_aij {
		row_rad_max_aij(double rad, double max_aij) : rad(rad), max_aij(max_aij) { }
		double rad;
		double max_aij;
	};

	const row_info compute_row_info(const affine& x, const double value) const;

	const row_rad_max_aij get_row_rad_max_aij(const affine& x) const;

	void set_col_bounds();

	void reset_col_arrays(int size);

	int  col_size() const;

	lp_impl* lp;
	int N_VARS;
	const double TINY;
	std::vector<std::vector<int> > pruning_indices;

	std::vector<affine>* v;

	std::vector<int>    col_index;
	std::vector<double> col_coeff;
};

}

#endif // LP_SOLVER_HPP_
