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

#include <stdint.h>
#include "lp_impl.hpp"
#include "glpk.h"

namespace asol {

class glpk_impl : public lp_impl {

public:

	glpk_impl();

private:

	virtual ~glpk_impl();

	virtual void reset();

	virtual void add_cols(int n);

	// index[1] ... index[length]
	virtual void add_eq_row(const int index[], const double value[], int length, double lb, double ub);

	virtual void set_col_bounds(int index, double lb, double ub);

	virtual void run_simplex();

	virtual void tighten_col_lb(int i, double& lb);

	virtual void tighten_col_ub(int i, double& ub);

	virtual int num_cols() const;

	virtual int num_rows() const;

	virtual col_status col_stat(int i) const;

	virtual double col_val(int i) const;

	virtual double col_lb(int i) const;

	virtual double col_ub(int i) const;

	virtual double col_dual_val(int i) const;

	virtual bool is_fixed(int index) const;

	virtual void dump(const char* file) const;

	virtual void show_iteration_count() const;

	static void free_environment();

	//===================================

	glpk_impl(const glpk_impl& );

	glpk_impl& operator=(const glpk_impl& );

	void init();

	void scale_prob();

	void warm_up_basis();

	void make_basis();

	void make_dual_feasible_basis();

	void set_col_dual_status(const int j);

	double solve_for(int index, int direction);

	//===================================

	glp_prob* lp;

	glp_smcp* parm;

	uint64_t previous_itr_count;
};

}

#endif // GLPK_IMPL_HPP_
