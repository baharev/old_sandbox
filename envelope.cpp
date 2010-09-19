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

#include <assert.h>
#include <ostream>
#include "envelope.hpp"
#include "lp_impl.hpp"
#include "exceptions.hpp"

lp_solver::lp_impl* asol::var::lp(new lp_solver::lp_impl());

namespace {

// Based on the C-XSC source code
void mult(const double xl, const double xu, const double yl, const double yu, double& zl, double& zu) {

	if (xl >=0) {                          /*  0 <= [x]                 */

		if (yl >=0)                        /*  0 <= [y]                 */
			zl=xl*yl;
		else                               /*  [y] <= 0  or  0 \in [y]  */
			zl=xu*yl;

		if (yu <=0)                        /*  [y] <= 0                 */
			zu=xl*yu;
		else                               /*  0 <= [y]  or  0 \in [y]  */
			zu=xu*yu;

	} else if (xu<=0) {                    /*  [x] <= 0                 */

		if (yu<=0)                         /*  [y] <= 0                 */
			zl=xu*yu;
		else                               /*  0 <= [y]  or  0 \in [y]  */
			zl=xl*yu;

		if (yl>=0)                         /*  0 <= [y]                 */
			zu=xu*yl;
		else                               /*  [y] <= 0  or  0 \in [y]  */
			zu=xl*yl;

	} else {                               /*  0 \in [x]                */

		if (yl>=0) {                       /*  0 <= [y]                 */
			zl=xl*yu;
			zu=xu*yu;
		} else if (yu<=0) {                /*  [y] <= 0                 */
			zl=xu*yl;
			zu=xl*yl;
		} else {                           /*  0 \in [x], 0 \in [y]     */
			const double lu = xl*yu;
			const double ul = xu*yl;
			const double ll = xl*yl;
			const double uu = xu*yu;
			zl=(lu<ul)?lu:ul;
			zu=(ll>uu)?ll:uu;
		}

	}

	assert(zl<=zu);
}

}

namespace asol {

void dbg_consistency(const var& a) {
	assert(a.lb <= a.ub);
	assert(a.index >= 1);
}

void dbg_consistency(const var& x, const var& y) {
	dbg_consistency(x);
	dbg_consistency(y);
}

void var::dump_lp(const char* file) {

	lp->dump(file);
}

var::var(double lb, double ub) : index(-1), lb(lb), ub(ub) {

	assert(lb <= ub);

	index = lp->add_col(lb, ub);
}

void var::fix_at(double val) {

	dbg_consistency(*this);

	if ((lb>val) || (val>ub)) {

		throw infeasible_problem();
	}

	lb = ub = val;

	lp->fix_col(index, val);
}

const var operator+(const var& x, const var& y) {

	dbg_consistency(x, y);

	double lb = x.lb + y.lb;
	double ub = x.ub + y.ub;

	var z(lb, ub);

	// x + y - z = 0
	var::lp->add_row(x.index, y.index, z.index);

	var::lp->tighten_col_bnds(z.index, z.lb, z.ub);

	return z;
}

const var operator*(const var& x, const var& y) {

	dbg_consistency(x, y);

	double lb =  1.0;
	double ub = -1.0;

	mult(x.lb, x.ub, y.lb, y.ub, lb, ub);

	var z(lb, ub);

	return z;
}

std::ostream& operator<<(std::ostream& os, const var& v) {

	return os << "[ " << v.lb << ", " << v.ub << "]" << std::flush;
}

}

