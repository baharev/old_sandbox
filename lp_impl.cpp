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

#include "lp_impl.hpp"
#include <iostream>

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
}

lp_impl::~lp_impl() {

	delete parm;

	glp_delete_prob(lp);

	glp_free_env();
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

bool lp_impl::simplex() {

	int ret = glp_simplex(lp, parm);

	if (ret) {
		clog << "Numerical problems in the simplex solver, code: " << ret << endl;
		return true;
	}
	return false;
}

void lp_impl::dump(const char* file) {

	glp_write_lp(lp, NULL, file);

}

}
