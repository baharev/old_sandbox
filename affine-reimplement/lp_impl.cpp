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

#include <iostream>
#include "lp_impl.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"

namespace asol {

lp_impl::lp_impl() {

	lp = glp_create_prob();

	parm = new glp_smcp;

	glp_init_smcp(parm);

	previous_itr_count = 0;

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

	parm->msg_lev = GLP_MSG_ALL;

	//parm->meth = GLP_DUAL;
}

void lp_impl::reset() {

	previous_itr_count += lpx_get_int_parm(lp, LPX_K_ITCNT);

	glp_erase_prob(lp);

	init();
}

void lp_impl::add_cols(int n) {

	glp_add_cols(lp, n);
}

void lp_impl::add_eq_row(const int index[], const double value[], int length, double lb, double ub) {

	using namespace std;

	//for (int i=1; i<=length; ++i) {

	//	cout << "index: " << index[i] << ", value: " << value[i] << endl;
	//}

	const int row_index = glp_add_rows(lp, 1);

	glp_set_mat_row(lp, row_index, length, index, value);

	const int row_type = (lb==ub) ? GLP_FX : GLP_DB;

	glp_set_row_bnds(lp, row_index, row_type, lb, ub);
}

void lp_impl::set_col_bounds(int index, const double lb, const double ub) {

	ASSERT(lb < ub);

	// FIXME Check if col bnds are better! Difficulty: lb=ub=0 by default

	glp_set_col_bnds(lp, index, GLP_DB, lb, ub);
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

void lp_impl::warm_up_basis() {

	int status = glp_warm_up(lp);

	if (status==0) {

		; // OK
	}
	else if (status==GLP_EBADB) {

		glp_std_basis(lp);

		std::cout << "Warning: bad basis!" << std::endl;

		warm_up_basis();
	}
	else if (status==GLP_ESING || status==GLP_ECOND) {

		throw numerical_problems();
	}
	else {

		ASSERT2(false,"Unknown status returned: "<<status);
	}
}

void lp_impl::make_basis() {

	warm_up_basis();

	if (parm->meth == GLP_PRIMAL) {

		; // Intentionally nothing at the moment
	}
	else {

		make_dual_feasible_basis();
	}
}

void lp_impl::make_dual_feasible_basis() {

	int status = glp_get_dual_stat(lp);

	if (status==GLP_FEAS || status==GLP_OPT) {

		return;
	}

	const int n = glp_get_num_cols(lp);

	for (int j=1; j<=n; ++j) {

		set_col_dual_status(j);
	}

	warm_up_basis();

	// TODO What happens if DUALP is used and falls back to std_basis and PRIMAL?
	status = glp_get_dual_stat(lp);

	ASSERT2(status==GLP_FEAS || status==GLP_OPT,"status: " << status);
}

void lp_impl::set_col_dual_status(const int j) {

	const int type = glp_get_col_type(lp, j);

	if (type!=GLP_FX) {

		const double d = glp_get_col_dual(lp, j);

		const int stat = (d>=0)?GLP_NL:GLP_NU;

		glp_set_col_stat(lp, j, stat);
	}
}

void lp_impl::run_simplex() {

	scale_prob();

	//make_basis(); // Redundant for PRIMAL, DUAL is unclear

	const int error_code = glp_simplex(lp, parm);

	if (error_code) {

		throw numerical_problems();
	}

	const int status = glp_get_status(lp);

	if (status==GLP_OPT) {

		;
	}
	else if (status==GLP_INFEAS || status==GLP_NOFEAS) { // TODO What is NOFEAS?

		throw infeasible_problem();
	}
	else {

		ASSERT(false); // Not clear how we could get kere
	}
}

double lp_impl::solve_for(int index, int direction) {

	glp_set_obj_dir(lp, direction);

	glp_set_obj_coef(lp, index, 1.0);

	if (parm->msg_lev >= GLP_MSG_ON) {

		std::cout << (direction==GLP_MIN?"MIN":"MAX") << std::endl;
	}

	try {

		run_simplex();
	}
	catch (...) {

		glp_set_obj_coef(lp, index, 0.0);

		throw;
	}

	return glp_get_col_prim(lp, index);
}

void lp_impl::tighten_col_lb(int i, double& lb) {

	ASSERT2(!is_fixed(i),"i: " << i);

	const double inf = solve_for(i, GLP_MIN);

	if (inf > lb) {

		lb = inf;
	}
}

void lp_impl::tighten_col_ub(int i, double& ub) {

	ASSERT2(!is_fixed(i),"i: " << i);

	const double sup = solve_for(i, GLP_MAX);

	if (sup < ub) {

		ub = sup;
	}
}

void lp_impl::dump(const char* file) const {

	glp_write_lp(lp, NULL, file);
}

void lp_impl::show_iteration_count() const {

	uint64_t count = previous_itr_count + lpx_get_int_parm(lp, LPX_K_ITCNT);

	std::cout << "Simplex iterations: " << count << std::endl;
}

int lp_impl::num_cols() const {

	return glp_get_num_cols(lp);
}

int lp_impl::num_rows() const {

	return glp_get_num_rows(lp);
}

col_status lp_impl::col_stat(int i) const {

	const int status = glp_get_col_stat(lp, i);

	col_status res;

	if (status==GLP_BS) {

		res = BASIC;
	}
	else if (status==GLP_NL) {

		res = NONBASIC_LB;
	}
	else if (status==GLP_NU) {

		res = NONBASIC_UB;
	}
	else {

		ASSERT2(false,"status: "<<status);
	}

	return res;
}

double lp_impl::col_val(int i) const {

	double val = glp_get_col_prim(lp, i);

	const double lb = glp_get_col_lb(lp, i);

	const double ub = glp_get_col_ub(lp, i);

	if (val < lb) {

		ASSERT(val+1.0e-4 > lb);

		val = lb;
	}
	else if (val > ub) {

		ASSERT(val-1.0e-4 < ub);

		val = ub;
	}

	return val;
}

double lp_impl::col_lb(int i) const {

	return glp_get_col_lb(lp, i);
}

double lp_impl::col_ub(int i) const {

	return glp_get_col_ub(lp, i);
}

bool lp_impl::is_fixed(int index) const {

	return glp_get_col_type(lp, index)==GLP_FX;
}

}
