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

// See the excellent documentation of GLPK for further details on the
// individual glp_ functions below

namespace lp_solver {

class lp_impl {

	public:

		lp_impl();

		int add_col(double lb, double ub);

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

		void remove_envelope(int index[5]);

		void fix_col(int index, double value);

		void set_bounds(int index, double lb, double ub);

		bool tighten_col_bnds(int i, double& lb, double& ub);

		bool col_type_db_or_fx(int index) const;

		void dump(const char* file) const;

		//int warm_up() { return glp_warm_up(lp); }

		double get_col_prim(int j) { return glp_get_col_prim(lp, j); }

		double get_obj_val()       { return glp_get_obj_val(lp); }

		int    get_it_cnt()        { return lpx_get_int_parm(lp, LPX_K_ITCNT); }

		void scale_prob(int flag) { glp_scale_prob(lp, flag); }

		void set_col_bnds(int j, int type, double lb, double ub) {

			glp_set_col_bnds(lp, j, type, lb, ub);
		}

		void set_row_bnds(int i, int type, double lb, double ub) {

			glp_set_row_bnds(lp, i, type, lb, ub);
		}

		void set_obj_coef(int j, double val) { glp_set_obj_coef(lp, j, val); }

		void set_msg_lev(const int level)    { parm->msg_lev = level; }

		void set_max_restart(const int);

		void std_basis()                     { glp_std_basis(lp); }

		int  term_out(int flag)              { return glp_term_out(flag); }

		~lp_impl();

	private:

		lp_impl(const lp_impl& );

		lp_impl& operator=(const lp_impl& );

		void throw_if_numerical_problems(int error, int line);

		void throw_if_infeasible(int status, int line);

		void throw_if_inconsistent_bnds(double lb, double ub, int line);

		void assert_col_type(int j, int line);

		void assert_value_within_bnds(int j, double value, int line);

		void assert_feasible_bounds(int j, double lb, double ub, int line);

		double solve_for(int index, int direction);

		void add_lu_row(double a, int x, int z, double c, int type);

		int add_lu_row(double a, int x, double b, int y, int z, double c, int type);

		void refresh(int index = 0);

		void reset_obj(int index);

		void write_back(int j, double inf, double sup, double& lb, double& ub);

		glp_prob* lp;

		glp_smcp* parm;
};

}

#endif /* LP_IMPL_HPP_ */
