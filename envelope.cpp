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

#include <ostream>
#include "envelope.hpp"
#include "lp_impl.hpp"

lp_solver::lp_impl* asol::var::lp(new lp_solver::lp_impl());

namespace asol {

class infeasible_problem {

};

class numerical_problems {

};

void var::dump_lp(const char* file) {

	lp->dump(file);
}

var::var(double lb, double ub) : index(-1), lb(lb), ub(ub) {

	index = lp->add_col(lb, ub);
}

void var::fix_at(double val) {

	if ((lb>val) || (val>ub)) {

		throw infeasible_problem();
	}

	lb = ub = val;

	lp->set_col_bnds(index, lp_solver::FX, lb, ub);

	//---

	bool error = lp->simplex();

	if (error) {

		throw numerical_problems();
	}

	int status = lp->get_status();

	if (status == lp_solver::NOFEAS) {

		throw infeasible_problem();
	}
	else if (status == lp_solver::OPT) {
		;
	}
	else {
		; // FIXME
	}

	//---
}

const var operator+(const var& x, const var& y) {

	double lb = x.lb + y.lb;
	double ub = x.ub + y.ub;

	var z(lb, ub);

	// x + y - z = 0
	var::lp->add_row(x.index, y.index, z.index);

	// TODO Finish!

	return z;
}

std::ostream& operator<<(std::ostream& os, const var& v){

	return os << "[ " << v.lb << ", " << v.ub << "]" << std::flush;
}

}

