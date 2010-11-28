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

#include <cmath>
#include <ostream>
// FIXME Remove and introduce logging
#include <iostream>
#include <assert.h>
#include "constants.hpp"
#include "envelope.hpp"
#include "lp_pair.hpp"
#include "exceptions.hpp"

lp_solver::lp_pair* asol::var::lp(new lp_solver::lp_pair);

namespace asol {

void dbg_consistency(const var& x, const var& y) {
	x.check_consistency();
	y.check_consistency();
}

void var::dump_lp(const char* file) {

	lp->dump(file);
}

var::var(const interval& bounds) : index(-1), range(bounds) {

	index = lp->add_col_nonbasic(range.inf(), range.sup());
}

var::var(double lb, double ub) : index(-1), range(interval(lb, ub)) {

	index = lp->add_col_nonbasic(lb, ub);
}

void var::fix_at(double val) {

	check_consistency();

	if (!range.contains(val)) {

		throw infeasible_problem();
	}

	range = interval(val);

	lp->fix_col(index, val);
}

const interval var::lp_tighten_col(bool& improved) const {

	check_consistency();

	double lb = range.inf();

	double ub = range.sup();

	improved = lp->tighten_col(index, lb, ub);

	return interval(lb, ub);
}

bool var::tighten_bounds() {

	bool improved = false;

	range = lp_tighten_col(improved);

	check_consistency();

	return improved;
}

const interval var::compute_bounds() const {

	bool dummy = false;

	return lp_tighten_col(dummy);
}

const var operator+(const var& x, const var& y) {

	var z(x.compute_bounds() + y.compute_bounds());

	// x + y - z = 0
	var::lp->add_add_row(x.index, y.index, z.index);

	z.tighten_bounds();

	return z;
}

const var operator+(const var& x, double y) {

	var z(x.compute_bounds() + y);

	// x - z = -y
	var::lp->add_shift_row(x.index, z.index, y);

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

	var z(x.compute_bounds() - y.compute_bounds());

	// x - y - z = 0
	var::lp->add_sub_row(x.index, y.index, z.index);

	z.tighten_bounds();

	return z;
}

const var operator-(double x, const var& y) {

	var z(x-y.compute_bounds());

	// z = x - y -> y + z = x
	var::lp->add_sum_row(y.index, z.index, x);

	z.tighten_bounds();

	return z;
}

bool var::contains_zero() const {

	check_consistency();

	return range.contains(0.0);
}

const var sqr(const var& x) {

	const interval x_range = x.compute_bounds();

	var z(sqr(x_range));

	const double xL = x_range.inf();
	const double xU = x_range.sup();

	// xL*xU <= (xL+xU)*x - z
	var::lp->add_lo_row(xL+xU, x.index, z.index, xL*xU);

	// 2*xL*x - z <= xL*xL
	var::lp->add_up_row(2*xL, x.index, z.index, xL*xL);

	// 2*xU*x - z <= xU*xU
	var::lp->add_up_row(2*xU, x.index, z.index, xU*xU);

	z.tighten_bounds();

	return z;
}

double y_eq(double x) {

	const double alpha = 3.55;

	return (alpha*x)/(1.0+(alpha-1.0)*x);
}

double y_derivative(double x) {

	const double alpha = 3.55;

	return alpha/std::pow(x*(alpha-1.0)+1.0, 2);
}

const var y_eq(const var& x) {

	const interval x_range = x.compute_bounds();

	const double xL = x_range.inf();
	const double xU = x_range.sup();

	const double yL = y_eq(xL);
	const double yU = y_eq(xU);

	var y(yL, yU);

	// m*x0-y0<=m*x-y
	double mL = y_derivative(xL);
	var::lp->add_lo_row(mL, x.index, y.index, mL*xL-yL);

	double mU = y_derivative(xU);
	var::lp->add_lo_row(mU, x.index, y.index, mU*xU-yU);

	double x_diam = x_range.diameter();
	double s = (x_diam>TOL_RANGE)?(yU-yL)/x_diam:y_derivative(x_range.midpoint());

	// s*xU-yU >= s*x-y
	var::lp->add_up_row(s, x.index, y.index, s*xU-yU);

	y.tighten_bounds();

	return y;
}

double H_liq(double x) {

	return 0.1667*std::exp(-1.087*x);
}

double H_liq_derivative(double x) {

	return -1.087*H_liq(x);
}

const var H_Liq(const var& x) {

	const interval x_range = x.compute_bounds();

	const double xL = x_range.inf();
	const double xU = x_range.sup();

	const double zL = H_liq(xU);
	const double zU = H_liq(xL);

	var z(zL, zU);

	// m*x0-y0>=m*x-y
	double mL = H_liq_derivative(xL);
	var::lp->add_up_row(mL, x.index, z.index, mL*xL-zU);

	double mU = H_liq_derivative(xU);
	var::lp->add_up_row(mU, x.index, z.index, mU*xU-zL);

	double x_diam = x_range.diameter();
	double s = (x_diam>TOL_RANGE)?(zL-zU)/x_diam:H_liq_derivative(x_range.midpoint());

	// s*xU-yU <= s*x-y
	var::lp->add_lo_row(s, x.index, z.index, s*xU-zL);

	z.tighten_bounds();

	return z;
}

double H_vap(double x) {

	return 0.1349*std::exp(-3.98*x) + 0.4397*std::exp(-0.088*x);
}

double H_vap_derivative(double x) {

	return (-3.98)*0.1349*std::exp(-3.98*x)+(-0.088)*0.4397*std::exp(-0.088*x);
}

const var H_Vap(const var& x) {

	const interval x_range = x.compute_bounds();

	const double xL = x_range.inf();
	const double xU = x_range.sup();

	const double zL = H_vap(xU);
	const double zU = H_vap(xL);

	var z(zL, zU);

	// m*x0-y0>=m*x-y
	double mL = H_vap_derivative(xL);
	var::lp->add_up_row(mL, x.index, z.index, mL*xL-zU);

	double mU = H_vap_derivative(xU);
	var::lp->add_up_row(mU, x.index, z.index, mU*xU-zL);

	double x_diam = x_range.diameter();
	double s = (x_diam>TOL_RANGE)?(zL-zU)/x_diam:H_vap_derivative(x_range.midpoint());

	// s*xU-yU <= s*x-y
	var::lp->add_lo_row(s, x.index, z.index, s*xU-zL);

	z.tighten_bounds();

	return z;
}

const var operator*(const var& x, const var& y) {

	dbg_consistency(x, y);

	if (x.index==y.index) {
		// FIXME You really should not write things like x*x ...
		return sqr(x);
	}

	const interval X = x.compute_bounds();
	const interval Y = y.compute_bounds();

	var z(X*Y);

	const double xL = X.inf();
	const double xU = X.sup();

	const double yL = Y.inf();
	const double yU = Y.sup();

	// TODO Make reset = false in lp_pair
	// TODO Use interval in lp_pair
	var::lp->add_mult_envelope(x.index, xL, xU, y.index, yL, yU, z.index, false);

	z.tighten_bounds();

	return z;
}

const var operator*(const double c, const var& x) {

	const interval X = x.compute_bounds();

	var z(c*X);

	var::lp->add_cx_row(c, x.index, z.index);

	z.tighten_bounds();

	return z;
}

const var operator*(const var& x, const double c) {

	return c*x;
}

void var::intersect(double lb, double ub) {

	tighten_bounds();

	bool improved = range.intersect(lb, ub);

	if (improved) {

		var::lp->set_bounds(index, range.inf(), range.sup());
	}
}

const var operator/(const var& x, const var& y) {

	assert(!y.contains_zero());

	if (x.index==y.index) {
		// FIXME You really should not write things like x/x ...
		return var(1,1);
	}

	interval Y = y.compute_bounds();
	interval X = x.compute_bounds();

	var z(X/Y);

	//z.intersect(lb, -12.0);
	//z.intersect(5.189, ub);
	//z.intersect(lb, -2.956);
	//z.intersect(lb, -1.0);

	bool improved = false;

	int counter = 0;

	interval Z = z.range;

	do {

		double diam_prev = Z.diameter();

		var::lp->add_mult_envelope(y.index, Y.inf(), Y.sup(), z.index, Z.inf(), Z.sup(), x.index, improved);

		improved =
		// FIXME Should check if z.ub <= col_ub in LP
		//var::lp_max->tighten_col_ub(z.index, z.ub);
		z.tighten_bounds();

		using namespace std;

		cout << endl << "z: " << z << ", pass: " << ++counter << endl;

		Z = z.range;

		if (((diam_prev-Z.diameter()))<1.0e-6*diam_prev)
			break;

		propagate_mult(X, Y, Z);
		// TODO Write back to lp
		//var::lp->set_bounds(x.index, X.inf(), X.sup());
		//var::lp->set_bounds(y.index, Y.inf(), Y.sup());

	} while (improved);

	return z;
}

std::ostream& operator<<(std::ostream& os, const var& v) {

	return os << v.range << std::flush;
}


void var::copy_bounds(double& lo, double& up) const {

	check_consistency();

	lo = range.inf();
	up = range.sup();
}

void var::reset() {

	lp->reset();
}

void var::release_all() {

	delete lp;
	lp = 0;
	lp_solver::lp_pair::free_environment();
}

void var::check_consistency() const {
	assert(range.inf() <= range.sup());
	assert(index >= 1);
	assert(var::lp->col_type_db_or_fx(index));
}

}
