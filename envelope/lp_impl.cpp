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
#include <cmath>
#include <assert.h>
#include "constants.hpp"
#include "exceptions.hpp"
#include "lp_impl.hpp"

using std::clog;
using std::cout;
using std::endl;
using std::fabs;

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

	init();
}

lp_impl::~lp_impl() {

	delete parm;

	glp_delete_prob(lp);
}

void lp_impl::free_environment() {

	glp_free_env();
}

void lp_impl::init() {

	glp_set_obj_dir(lp, GLP_MIN);

	parm->presolve = GLP_OFF;

	parm->msg_lev = GLP_MSG_OFF;

	//parm->meth = GLP_DUAL;

	dual_feasible = true;

	assert((dual_feasible==true && parm->meth==GLP_PRIMAL)||
		  ((dual_feasible==false&& parm->meth==GLP_DUAL)) );
}

void lp_impl::reset() {

	glp_erase_prob(lp);

	init();
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

void lp_impl::assert_col_type(int j, int line) {

	int type = glp_get_col_type(lp, j);

	if (!col_type_db_or_fx(j)) {
		clog << "Error: col " << j << " is of type " << type << "; ";
		clog << __FILE__ ", line " << line << endl;
		throw asol::assertion_error();
	}
}

double AbsMax(double x, double y) {

	double X = fabs(x);

	double Y = fabs(y);

	return (X<Y)?Y:X;
}

void lp_impl::assert_value_within_bnds(int j, double value, int line) {

	assert_col_type(j, line);

	double lb = glp_get_col_lb(lp, j);
	double ub = glp_get_col_ub(lp, j);

	if (value < lb || value > ub) {
		clog << "Error: " << value << " is not in range ";
		clog << "[ " << lb << ", " << ub << "]" << endl;
		// This should have been checked in envelope.cpp
		throw asol::assertion_error();
	}
}

void lp_impl::assert_feasible_bounds(int j, double l, double u, int line) {

	assert_col_type(j, line);

	double lb = glp_get_col_lb(lp, j);
	double ub = glp_get_col_ub(lp, j);

	if (lb > l || ub < u) {

		throw asol::assertion_error();
	}
}

void lp_impl::set_max_restart(const int limit) {

	set_restart_limit(limit);
}

int lp_impl::add_col_nonbasic_on_lb(double lb, double ub) {

	return add_new_col(lb, ub, GLP_NL);
}

int lp_impl::add_col_nonbasic_on_ub(double lb, double ub) {

	return add_new_col(lb, ub, GLP_NU);
}

int lp_impl::add_new_col(double lb, double ub, int stat) {

	int j = glp_add_cols(lp, 1);

	int type = (lb == ub)?GLP_FX:GLP_DB;

	glp_set_col_bnds(lp, j, type, lb, ub);

	glp_set_col_stat(lp, j, stat);

	return j;
}

// x + y = c
void lp_impl::add_sum_row(int x, int y, double c) {

	int i = add_new_row(GLP_FX, c);

	int ind[] = { 0, x, y };

	double val[] = { 0.0, 1.0, 1.0 };

	glp_set_mat_row(lp, i, 2, ind, val);
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
}

void lp_impl::get_row_status(const int rows[5], int stat[5]) const {

	for (int i=1; i<=4; ++i) {
		stat[i] = glp_get_row_stat(lp, rows[i]);
	}
}

void lp_impl::set_row_status(const int rows[5], const int stat[5]) {

	for (int i=1; i<=4; ++i) {
		glp_set_row_stat(lp, rows[i], stat[i]);
	}
}

int lp_impl::add_new_row(int type, double bound) {

	int i = glp_add_rows(lp, 1);

	glp_set_row_bnds(lp, i, type, bound, bound);

	glp_set_row_stat(lp, i, GLP_BS);

	return i;
}

// (c <=) ax - z (<= c)
void lp_impl::add_lu_row(double a, int x, int z, double c, int type) {

	int i = add_new_row(type, c);

	int ind[] = { 0, x, z };

	double val[] = { 0.0, a, -1.0 };

	glp_set_mat_row(lp, i, 2, ind, val);
}

// (c <=) ax + by - z (<= c)
int lp_impl::add_lu_row(double a, int x, double b, int y, int z, double c, int type) {

	int i = add_new_row(type, c);

	int ind[] = { 0, x, y, z };

	double val[] = { 0.0, a, b, -1.0 };

	glp_set_mat_row(lp, i, 3, ind, val);

	return i;
}

void lp_impl::fix_col(int index, double value) {

	assert_value_within_bnds(index, value, __LINE__);

	glp_set_col_bnds(lp, index, GLP_FX, value, value);

	refresh();
}

void lp_impl::bounds_to_be_set(int index, double& l, double& u) {

	assert(l<=u);

	const double lb = glp_get_col_lb(lp, index);
	const double ub = glp_get_col_ub(lp, index);

	assert(lb<=ub);

	if (l < lb) {
		l = lb;
	}

	if (ub < u) {
		u = ub;
	}

	if (l>u) {
		throw asol::infeasible_problem();
	}
}

void lp_impl::set_bounds(int index, double lb, double ub) {

	bounds_to_be_set(index, lb, ub);

	assert_feasible_bounds(index, lb, ub, __LINE__);

	int type = (lb < ub)?GLP_DB:GLP_FX;

	glp_set_col_bnds(lp, index, type, lb, ub);

	refresh();
}

void lp_impl::make_dual_feas_basis() {

	glp_std_basis(lp);

	glp_warm_up(lp);

	const int n = glp_get_num_cols(lp);

	for (int j=1; j<=n; ++j) {

		const double d = glp_get_col_dual(lp, j);

		const int stat = (d>=0)?GLP_NL:GLP_NU;

		glp_set_col_stat(lp, j, stat);
	}

	glp_warm_up(lp);

	dual_feasible = true;
}

void lp_impl::scale_prob() {

	if (parm->msg_lev < GLP_MSG_ON) {

		glp_term_out(GLP_OFF);
	}

	glp_scale_prob(lp, GLP_SF_EQ);

	if (parm->msg_lev < GLP_MSG_ON) {

		glp_term_out(GLP_ON);
	}
}

void lp_impl::refresh(int index) {

	if (!dual_feasible) {

		make_dual_feas_basis();
	}

	scale_prob();

	int error_code = glp_simplex(lp, parm);

	if (index!=0) {

		reset_obj(index); // TODO Not in solve for because of throw, necessary?
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

	if (parm->msg_lev >= GLP_MSG_ON) {

		cout << (direction==GLP_MIN?"MIN":"MAX") << endl;
	}

	refresh(index);

	return glp_get_col_prim(lp, index);
}

bool lp_impl::is_fixed(int index) {

	return glp_get_col_type(lp, index)==GLP_FX;
}

bool lp_impl::tighten_col_lb(int index, double& lb) {

	assert(!is_fixed(index));

	bool improved = false;

	const double inf = solve_for(index, GLP_MIN);

	if (inf>lb) {
		lb = inf;
		improved = true;
	}

	return improved;
}

bool lp_impl::tighten_col_ub(int index, double& ub) {

	assert(!is_fixed(index));

	bool improved = false;

	const double sup = solve_for(index, GLP_MAX);

	if (sup<ub) {
		ub = sup;
		improved = true;
	}

	return improved;
}

/*
bool lp_impl::tighten_col_bnds(int index, double& lb, double& ub) {

	const bool improved_lb = tighten_col_lb(index, lb);

	const bool improved_ub = tighten_col_ub(index, ub);

	throw_if_inconsistent_bnds(lb, ub, __LINE__);

	return improved_lb || improved_ub;
}
*/

void lp_impl::dump(const char* file) const {

	glp_write_lp(lp, NULL, file);

}

}
