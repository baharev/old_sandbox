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

	glp_erase_prob(lp);
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

void lp_impl::scale_prob() {

	if (parm->msg_lev < GLP_MSG_ON) {

		glp_term_out(GLP_OFF);
	}

	glp_scale_prob(lp, GLP_SF_EQ);

	if (parm->msg_lev < GLP_MSG_ON) {

		glp_term_out(GLP_ON);
	}
}

void lp_impl::check_feasibility() {

	scale_prob();

	glp_std_basis(lp);

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

}
