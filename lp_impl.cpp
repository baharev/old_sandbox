//=============================================================================
//
//  This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
//  Copyright (C) 2010  Ali Baharev
//  All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//=============================================================================

#include <iostream>
#include "exceptions.hpp"
#include "lp_impl.hpp"

using std::clog;
using std::endl;

//#define HACKED_GLPK
#ifdef HACKED_GLPK
extern "C" void set_restart_limit(const int limit);
#else
static void set_restart_limit(const int limit) { }
#endif

namespace lp_solver {

lp_impl::lp_impl() {

	lp = glp_create_prob();

	parm = new glp_smcp;

	glp_init_smcp(parm);

	glp_set_obj_dir(lp, GLP_MIN);
}

lp_impl::~lp_impl() {

	delete parm;

	glp_delete_prob(lp);

	glp_free_env();
}

void lp_impl::throw_if_numerical_problems(int error, int line) {

	if (error != 0) {

		clog << "Numerical problems, code: " << error << "; ";
		clog << __FILE__ ", line " << line << endl;

		throw asol::numerical_problems();
	}
}

void lp_impl::throw_if_infeasible(int status, int line) {

	if (status == GLP_OPT) {
		return;
	}

	if (status == GLP_NOFEAS) {
		throw asol::infeasible_problem();
	}

	clog << "Unexpected status: " << status << "; " __FILE__ ", ";
	clog << "line " << line << endl;

	throw asol::numerical_problems();
}

void lp_impl::throw_if_inconsistent_bnds(double lb, double ub, int line) {

	if (lb > ub) {

		clog << "Inconsistent bounds, lb > ub: " << lb << " > " << ub << "; ";
		clog << __FILE__ ", line " << line << endl;

		throw asol::numerical_problems();
	}
}

void lp_impl::set_max_restart(const int limit) {

	set_restart_limit(limit);
}

int lp_impl::add_col(double lb, double ub) {

	int j = glp_add_cols(lp, 1);

	glp_set_col_bnds(lp, j, GLP_DB, lb, ub);

	glp_set_col_stat(lp, j, GLP_NL);

	return j;
}

void lp_impl::add_row(int x, int y, int z) {

	int i = glp_add_rows(lp, 1);

	glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);

	int ind[] = { 0, x, y, z };

	double val[] = { 0.0, 1.0, 1.0, -1.0 };

	glp_set_mat_row(lp, i, 3, ind, val);

	glp_set_row_stat(lp, i, GLP_BS);
}

void lp_impl::fix_col(int index, double value) {

	glp_set_col_bnds(lp, index, GLP_FX, value, value);

	refresh();
}

void lp_impl::refresh(int index) {

	int error_code = glp_simplex(lp, parm);

	if (index!=0) {

		reset_obj(index);
	}

	throw_if_numerical_problems(error_code, __LINE__);

	int status = glp_get_status(lp);

	throw_if_infeasible(status, __LINE__);
}

void lp_impl::reset_obj(int index) {

	glp_set_obj_dir(lp, GLP_MIN);

	glp_set_obj_coef(lp, index, 0.0);
}

void lp_impl::write_back(int j, double inf, double sup, double& lb, double& ub) {

	bool improved = false;

	if (inf>lb) {
		lb = inf;
		improved = true;
	}

	if (sup<ub) {
		ub = sup;
		improved = true;
	}

	if (improved) {
		glp_set_col_bnds(lp, j, GLP_DB, lb, ub);
	}
}

double lp_impl::solve_for(int index, int direction) {

	glp_set_obj_dir(lp, direction);

	glp_set_obj_coef(lp, index, 1.0);

	refresh(index);

	return glp_get_col_prim(lp, index);
}

void lp_impl::tighten_col_bnds(int index, double& lb, double& ub) {

	const double inf = solve_for(index, GLP_MIN);

	const double sup = solve_for(index, GLP_MAX);

	throw_if_inconsistent_bnds(inf, sup, __LINE__);

	write_back(index, inf, sup, lb, ub);
}

void lp_impl::dump(const char* file) {

	glp_write_lp(lp, NULL, file);

}

}
