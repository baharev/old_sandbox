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

#ifndef LP_IMPL_HPP_
#define LP_IMPL_HPP_

namespace asol {

enum col_status {
	BASIC,
	NONBASIC_LB,
	NONBASIC_UB
};

class lp_impl {

public:

	virtual void reset() = 0;

	virtual void add_cols(int n) = 0; // TODO Merge with reset?

	// index[1] ... index[length]
	virtual void add_eq_row(const int index[], const double value[], int length, double lb, double ub) = 0;

	virtual void set_col_bounds(int index, double lb, double ub) = 0;

	virtual void run_simplex() = 0;

	virtual void tighten_col_lb(int i, double& lb) = 0; // TODO Return the new bound instead?

	virtual void tighten_col_ub(int i, double& ub) = 0;

	virtual int num_cols() const = 0;

	virtual int num_rows() const = 0;

	virtual col_status col_stat(int i) const = 0;

	virtual double col_val(int i) const = 0;

	virtual double col_lb(int i) const = 0;

	virtual double col_ub(int i) const = 0;

	virtual double col_dual_val(int i) const = 0;

	virtual bool is_fixed(int index) const = 0;

	virtual void dump(const char* file) const = 0;

	virtual void show_iteration_count() const = 0;

	virtual ~lp_impl() = 0;

protected:

	lp_impl() { }
};

inline lp_impl::~lp_impl() { }

}

#endif // LP_IMPL_HPP_
