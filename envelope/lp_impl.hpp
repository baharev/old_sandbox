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

#ifndef LP_IMPL_HPP_
#define LP_IMPL_HPP_

#include "glpk.h"

namespace asol {

class lp_impl {

	public:

		lp_impl();

		void reset();

		int add_col_nonbasic_on_lb(double lb, double ub);

		int add_col_nonbasic_on_ub(double lb, double ub);

		// x + y = c
		void add_sum_row(int x, int y, double c);

		// x - z = -y
		void add_shift_row(int x, int z, double y);

		// c*x - z = 0
		void add_cx_row(double c, int x, int z);

		// x + y - z = 0
		void add_add_row(int x, int y, int z);

		// x - y - z = 0
		void add_sub_row(int x, int y, int z);

		// c <= ax - z
		void add_lo_row(double a, int x, int z, double c);

		// ax - z <= c
		void add_up_row(double a, int x, int z, double c);

		// c <= ax + by - z
		int add_lo_row(double a, int x, double b, int y, int z, double c);

		// ax + by - z <= c
		int add_up_row(double a, int x, double b, int y, int z, double c);

		// sum c*x = 0
		void add_lin_con(double val, const double c[], const int x[], int size);

		void remove_envelope(int index[5]);

		void get_row_status(const int rows[5], int stat[5]) const;

		void set_row_status(const int rows[5], const int stat[5]);

		void fix_col(int index, double value);

		bool is_fixed(int index);

		void set_bounds(int index, double lb, double ub);

		void refresh_basis();

		double get_col_val(int index);

		bool tighten_col_lb(int i, double& lb);

		bool tighten_col_ub(int i, double& ub);

		int n_cols();

		double col_lb(int i);

		double col_ub(int i);

		bool col_type_db_or_fx(int index) const;

		void dump(const char* file) const;

		void set_max_restart(const int);

		~lp_impl();

		static void free_environment();

	private:

		lp_impl(const lp_impl& );

		lp_impl& operator=(const lp_impl& );

		void init();

		void throw_if_numerical_problems(int error, int line);

		void throw_if_infeasible(int status, int line);

		void throw_if_inconsistent_bnds(double lb, double ub, int line);

		void assert_col_type(int j, int line);

		void assert_value_within_bnds(int j, double value, int line);

		void assert_feasible_bounds(int j, double lb, double ub, int line);

		double solve_for(int index, int direction);

		int add_new_col(double lb, double ub, int status);

		int add_new_row(int type, double bound);

		void add_lu_row(double a, int x, int z, double c, int type);

		int add_lu_row(double a, int x, double b, int y, int z, double c, int type);

		void make_dual_feas_basis();

		void scale_prob();

		void refresh(int index = 0);

		void reset_obj(int index);

		void bounds_to_be_set(int index, double& l, double& u);

		glp_prob* lp;

		glp_smcp* parm;
		 // TODO Find a better name, in case of primal, it is misleading
		bool dual_feasible;
};

bool too_narrow(double lb, double ub);

}

#endif /* LP_IMPL_HPP_ */
