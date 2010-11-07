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
#include <string>
// FIXME Remove and introduce logging
#include <iostream>
#include "envelope.hpp"
#include "lp_impl.hpp"
#include "exceptions.hpp"


lp_solver::lp_impl* asol::var::lp_min(new lp_solver::lp_impl());
lp_solver::lp_impl* asol::var::lp_max(new lp_solver::lp_impl());

namespace {

// Based on the C-XSC source code
void mult(const double xl, const double xu, const double yl, const double yu, double& zl, double& zu) {

	assert(xl<=xu);
	assert(yl<=yu);

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

// Based on the C-XSC source code
void division(const double xl, const double xu, const double yl, const double yu, double& zl, double& zu) {

	assert(xl<=xu);
	assert(yl<=yu);

	if (yl>0) {

		if (xl>=0)
			zl=xl/yu;
		else
			zl=xl/yl;

		if (xu<=0)
			zu=xu/yu;
		else
			zu=xu/yl;

	}
	else if (yu<0) {

		if (xu<=0)
			zl=xu/yl;
		else
			zl=xu/yu;

		if (xl>=0)
			zu=xl/yl;
		else
			zu=xl/yu;
	}
	else {
		assert(false);
	}

	assert(zl<=zu);
}

void sort(double& lb, double& ub) {

	if (lb > ub) {
		double temp = lb;
		lb = ub;
		ub = temp;
	}
}

}

namespace asol {

void dbg_consistency(const var& x, const var& y) {
	x.check_consistency();
	y.check_consistency();
}

void var::dump_lp(const char* file) {

	using std::string;

	string fmin(file);
	string fmax(file);

	fmin += "_min";
	fmax += "_max";

	lp_min->dump(fmin.c_str());
	lp_max->dump(fmax.c_str());
}

var::var(double lb, double ub) : index(-1), lb(lb), ub(ub) {

	assert(lb <= ub);

	int index_min = lp_min->add_col_nonbasic_on_lb(lb, ub);
	int index_max = lp_max->add_col_nonbasic_on_ub(lb, ub);
	assert(index_min == index_max);
	index = index_min;
}

void var::fix_at(double val) {

	check_consistency();

	if ((lb>val) || (val>ub)) {

		throw infeasible_problem();
	}

	lb = ub = val;

	lp_min->fix_col(index, val);
	lp_max->fix_col(index, val);
}

bool var::tighten_bounds() {

	check_consistency();

	bool min_improved = lp_min->tighten_col_lb(index, lb);
	bool max_improved = lp_max->tighten_col_ub(index, ub);

	return min_improved || max_improved;
}

const var operator+(const var& x, const var& y) {

	dbg_consistency(x, y);

	double lb = x.lb + y.lb;
	double ub = x.ub + y.ub;

	var z(lb, ub);

	// x + y - z = 0
	var::lp_min->add_add_row(x.index, y.index, z.index);
	var::lp_max->add_add_row(x.index, y.index, z.index);

	z.tighten_bounds();

	return z;
}

const var operator+(const var& x, double y) {

	x.check_consistency();

	var z(x.lb+y, x.ub+y);

	// x - z = -y
	var::lp_min->add_shift_row(x.index, z.index, y);
	var::lp_max->add_shift_row(x.index, z.index, y);

	z.tighten_bounds();

	return z;
}

const var operator+(double x, const var& y) {

	return y+x;
}

const var operator-(const var& x, double y) {

	return x+(-y);
}

const var operator-(const var& x, const var& y) {

	dbg_consistency(x, y);

	double lb = x.lb - y.ub;
	double ub = x.ub - y.lb;

	var z(lb, ub);

	// x - y - z = 0
	var::lp_min->add_sub_row(x.index, y.index, z.index);
	var::lp_max->add_sub_row(x.index, y.index, z.index);

	z.tighten_bounds();

	return z;
}

bool var::contains_zero() const {

	check_consistency();

	return (lb<=0)&&(0<=ub);
}

void var::lp_add_lo_row(double a, int x, int z, double c) {

	lp_min->add_lo_row(a, x, z, c);
	lp_max->add_lo_row(a, x, z, c);
}

void var::lp_add_up_row(double a, int x, int z, double c) {

	lp_min->add_up_row(a, x, z, c);
	lp_max->add_up_row(a, x, z, c);
}

int var::lp_add_lo_row(double a, int x, double b, int y, int z, double c) {

	int i_min = lp_min->add_lo_row(a, x, b, y, z, c);
	int i_max = lp_max->add_lo_row(a, x, b, y, z, c);
	assert(i_min == i_max);
	return i_min;
}

int var::lp_add_up_row(double a, int x, double b, int y, int z, double c) {

	int i_min = lp_min->add_up_row(a, x, b, y, z, c);
	int i_max = lp_max->add_up_row(a, x, b, y, z, c);
	assert(i_min == i_max);
	return i_min;
}

const var sqr(const var& x) {

	x.check_consistency();

	double lb = x.lb*x.lb;
	double ub = x.ub*x.ub;

	sort(lb, ub);

	if (x.contains_zero()) {

		lb = 0.0;
	}

	var z(lb, ub);

	// xL*xU <= (xL+xU)*x - z
	var::lp_add_lo_row(x.lb+x.ub, x.index, z.index, x.lb*x.ub);

	// 2*xL*x - z <= xL*xL
	var::lp_add_up_row(2*x.lb, x.index, z.index, x.lb*x.lb);

	// 2*xU*x - z <= xU*xU
	var::lp_add_up_row(2*x.ub, x.index, z.index, x.ub*x.ub);

	z.tighten_bounds();

	return z;
}

void add_mult_envelope(const var& x, const var& y, const var& z, bool reset) {

	static int rows[] = { -1, -1, -1, -1, -1 };
	static int stat_min[] = { -1, -1, -1, -1, -1 };
	static int stat_max[] = { -1, -1, -1, -1, -1 };

	if (reset) {
		var::lp_min->get_row_status(rows, stat_min);
		var::lp_max->get_row_status(rows, stat_max);
		var::lp_min->remove_envelope(rows);
		var::lp_max->remove_envelope(rows);
	}

	// yL*xU <= yL*x + xU*y - z
	rows[1] = var::lp_add_lo_row(y.lb, x.index, x.ub, y.index, z.index, y.lb*x.ub);

	// yU*xL <= yU*x + xL*y - z
	rows[2] = var::lp_add_lo_row(y.ub, x.index, x.lb, y.index, z.index, y.ub*x.lb);

	// yL*x + xL*y - z <= yL*xL
	rows[3] = var::lp_add_up_row(y.lb, x.index, x.lb, y.index, z.index, y.lb*x.lb);

	// yU*x + xU*y - z <= yU*xU
	rows[4] = var::lp_add_up_row(y.ub, x.index, x.ub, y.index, z.index, y.ub*x.ub);

	if (reset) {
		var::lp_min->set_row_status(rows, stat_min);
		var::lp_max->set_row_status(rows, stat_max);
	}
}

const var operator*(const var& x, const var& y) {

	dbg_consistency(x, y);

	if (x.index==y.index) {
		// FIXME You really should not write things like x*x ...
		return sqr(x);
	}

	double lb =  1.0;
	double ub = -1.0;

	mult(x.lb, x.ub, y.lb, y.ub, lb, ub);

	var z(lb, ub);

	add_mult_envelope(x, y, z);

	z.tighten_bounds();

	return z;
}

const var operator*(const double c, const var& x) {

	x.check_consistency();

	double lb = c*x.lb;
	double ub = c*x.ub;

	sort(lb, ub);

	var z(lb, ub);

	var::lp_min->add_cx_row(c, x.index, z.index);
	var::lp_max->add_cx_row(c, x.index, z.index);

	z.tighten_bounds();

	return z;
}

const var operator*(const var& x, const double c) {

	return c*x;
}

void var::intersect(double l, double u) {

	check_consistency();

	bool improved = false;

	if (l>lb) {
		lb = l;
		improved = true;
	}

	if (u<ub) {
		ub = u;
		improved = true;
	}

	if (lb>ub)
		throw infeasible_problem();

	if (improved) {
		var::lp_min->set_bounds(index, lb, ub);
		var::lp_max->set_bounds(index, lb, ub);
	}
}

// z = x*y
void var::propagate(var& x, var& y) {

	check_consistency();
	dbg_consistency(x, y);

	if (!x.contains_zero()) { // y = z/x
		double y_lb =  1.0;
		double y_ub = -1.0;
		division(lb, ub, x.lb, x.ub, y_lb, y_ub);
		y.intersect(y_lb, y_ub);
	}

	if (!y.contains_zero()) {  // x = z/y
		double x_lb =  1.0;
		double x_ub = -1.0;
		division(lb, ub, y.lb, y.ub, x_lb, x_ub);
		x.intersect(x_lb, x_ub);
	}

	double z_lb =  1.0; // z = x*y
	double z_ub = -1.0;
	mult(x.lb, x.ub, y.lb, y.ub, z_lb, z_ub);
	intersect(z_lb, z_ub);
}

const var operator/(var& x, var& y) {

	dbg_consistency(x, y);

	assert(!y.contains_zero());

	if (x.index==y.index) {
		// FIXME You really should not write things like x/x ...
		return var(1,1);
	}

	double lb =  1.0;
	double ub = -1.0;

	division(x.lb, x.ub, y.lb, y.ub, lb, ub);

	var z(lb, ub);

	bool improved = false;

	do {

		add_mult_envelope(y, z, x, improved);

		improved =
		// FIXME Should check if z.ub <= col_ub in LP
		//var::lp_max->tighten_col_ub(z.index, z.ub);
		z.tighten_bounds();

		//z.intersect(5.189, 482);

		std::cout << std::endl << "z: " << z << std::endl;

		x.propagate(y, z);

	} while (improved);

	return z;
}

std::ostream& operator<<(std::ostream& os, const var& v) {

	return os << "[ " << v.lb << ", " << v.ub << "]" << std::flush;
}

void var::check_consistency() const {
	assert(lb <= ub);
	assert(index >= 1);
	assert(var::lp_min->col_type_db_or_fx(index));
	assert(var::lp_max->col_type_db_or_fx(index));
}

}

