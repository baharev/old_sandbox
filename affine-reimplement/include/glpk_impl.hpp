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

#ifndef GLPK_IMPL_HPP_
#define GLPK_IMPL_HPP_

#include "stdint.h"
#include "glpk.h"

namespace asol {

enum col_status {
	BASIC,
	NONBASIC_LB,
	NONBASIC_UB
};

class lp_impl {

public:

	lp_impl();

	~lp_impl();

	void reset();

	void add_cols(int n);

	// index[1] ... index[length]
	void add_eq_row(const int index[], const double value[], int length, double lb, double ub);

	void set_col_bounds(int index, double lb, double ub);

	void run_simplex();

	void tighten_col_lb(int i, double& lb);

	void tighten_col_ub(int i, double& ub);

	int num_cols() const;

	int num_rows() const;

	col_status col_stat(int i) const;

	double col_val(int i) const;

	double col_lb(int i) const;

	double col_ub(int i) const;

	double col_dual_val(int i) const;

	bool is_fixed(int index) const;

	void dump(const char* file) const;

	void show_iteration_count() const;

	static void free_environment();

private:

	lp_impl(const lp_impl& );

	lp_impl& operator=(const lp_impl& );

	void init();

	void scale_prob();

	void warm_up_basis();

	void make_basis();

	void make_dual_feasible_basis();

	void set_col_dual_status(const int j);

	double solve_for(int index, int direction);

	glp_prob* lp;

	glp_smcp* parm;

	uint64_t previous_itr_count;
};

}

#endif // GLPK_IMPL_HPP_
