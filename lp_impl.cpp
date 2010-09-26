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

	parm->presolve = GLP_OFF;

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

bool lp_impl::col_type_db_or_fx(int j) const {

	int type = glp_get_col_type(lp, j);

	return (type==GLP_DB) || (type==GLP_FX);
}

void lp_impl::assert_value_within_bnds(int j, double value, int line) {

	int type = glp_get_col_type(lp, j);

	if (!col_type_db_or_fx(j)) {
		clog << "Error: col " << j << " is of type " << type << "; ";
		clog << __FILE__ ", line " << line << endl;
		throw asol::assertion_error();
	}

	double lb = glp_get_col_lb(lp, j);
	double ub = glp_get_col_ub(lp, j);

	if (value < lb || value > ub) {
		// This should have been checked in envelope.cpp
		throw asol::assertion_error();
	}
}

void lp_impl::set_max_restart(const int limit) {

	set_restart_limit(limit);
}

int lp_impl::add_col(double lb, double ub) {

	int j = glp_add_cols(lp, 1);

	int type = (lb == ub)?GLP_FX:GLP_DB;

	glp_set_col_bnds(lp, j, type, lb, ub);

	glp_set_col_stat(lp, j, GLP_NL);

	return j;
}

// x - z = -y
void lp_impl::add_shift_row(int x, int z, double y) {

	add_lu_row(1.0, x, z, -y, GLP_FX);
}

// c*x - z = 0
void lp_impl::add_cx_row(double c, int x, int z) {

	add_lu_row(c, x, z, 0, GLP_FX);
}

// c <= ax - z
void lp_impl::add_lo_row(double a, int x, int z, double c) {

	add_lu_row(a, x, z, c, GLP_LO);
}

// ax - z <= c
void lp_impl::add_up_row(double a, int x, int z, double c) {

	add_lu_row(a, x, z, c, GLP_UP);
}

// c <= ax + by - z
int lp_impl::add_lo_row(double a, int x, double b, int y, int z, double c) {

	return add_lu_row(a, x, b, y, z, c, GLP_LO);
}

// ax + by - z <= c
int lp_impl::add_up_row(double a, int x, double b, int y, int z, double c) {

	return add_lu_row(a, x, b, y, z, c, GLP_UP);
}

// x + y - z = 0
void lp_impl::add_add_row(int x, int y, int z) {

	add_lu_row(1.0, x, 1.0, y, z, 0.0, GLP_FX);
}

// x - y - z = 0
void lp_impl::add_sub_row(int x, int y, int z) {

	add_lu_row(1.0, x, -1.0, y, z, 0.0, GLP_FX);
}

void lp_impl::remove_envelope(int index[5]) {

	glp_del_rows(lp, 4, index);
	// TODO Figure out a more efficient way...
	glp_std_basis(lp);
}

// (c <=) ax - z (<= c)
void lp_impl::add_lu_row(double a, int x, int z, double c, int type) {

	int i = glp_add_rows(lp, 1);

	glp_set_row_bnds(lp, i, type, c, c);

	int ind[] = { 0, x, z };

	double val[] = { 0.0, a, -1.0 };

	glp_set_mat_row(lp, i, 2, ind, val);

	glp_set_row_stat(lp, i, GLP_BS);
}

// (c <=) ax + by - z (<= c)
int lp_impl::add_lu_row(double a, int x, double b, int y, int z, double c, int type) {

	int i = glp_add_rows(lp, 1);

	glp_set_row_bnds(lp, i, type, c, c);

	int ind[] = { 0, x, y, z };

	double val[] = { 0.0, a, b, -1.0 };

	glp_set_mat_row(lp, i, 3, ind, val);

	glp_set_row_stat(lp, i, GLP_BS);

	return i;
}

void lp_impl::fix_col(int index, double value) {

	assert_value_within_bnds(index, value, __LINE__);

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

double lp_impl::solve_for(int index, int direction) {

	glp_set_obj_dir(lp, direction);

	glp_set_obj_coef(lp, index, 1.0);

	refresh(index);

	return glp_get_col_prim(lp, index);
}

bool lp_impl::tighten_col_bnds(int index, double& lb, double& ub) {

	const double inf = solve_for(index, GLP_MIN);

	const double sup = solve_for(index, GLP_MAX);

	throw_if_inconsistent_bnds(inf, sup, __LINE__);

	bool improved = false;

	if (inf>lb) {
		lb = inf;
		improved = true;
	}

	if (sup<ub) {
		ub = sup;
		improved = true;
	}

	return improved;
}

void lp_impl::dump(const char* file) const {

	glp_write_lp(lp, NULL, file);

}

}
