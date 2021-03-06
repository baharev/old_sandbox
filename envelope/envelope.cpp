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
#include "dag.hpp"
#include "lp_pair.hpp"
#include "exceptions.hpp"

// FIXME Clean up this mess!!!
#include "lp_pruning.hpp"

namespace asol {

lp_pair* var::lp(new lp_pair);
dag* var::ia_dag(new dag);

void dbg_consistency(const var& x, const var& y) {
	x.check_consistency();
	y.check_consistency();
}

void var::dump_lp(const char* file) {

	lp->dump(file);
}

var::var() : index(-1) { }

var::var(const interval& bounds) : index(-1) {

	index = lp->add_col_nonbasic(bounds.inf(), bounds.sup());

	ia_dag->add(index, bounds);
}

var::var(double lb, double ub) : index(-1) {

	index = lp->add_col_nonbasic(lb, ub);

	ia_dag->add(index, lb, ub);
}

bool var::intersect_in_dag(const interval& range) {

	return ia_dag->intersect(index, range);
}

const interval var::bounds() const {

	return ia_dag->bounds(index);
}

void var::fix_at(double val) {

	check_consistency();

	intersect_in_dag(interval(val));

	lp->fix_col(index, val);
}

const interval var::lp_tighten_col(bool& improved) const {

	check_consistency();

	interval range = bounds();

	double lb = range.inf();

	double ub = range.sup();

	improved = lp->tighten_col(index, lb, ub);

	return interval(lb, ub);
}

bool var::tighten_bounds() {

	bool improved = false;

	interval range = lp_tighten_col(improved);

	intersect_in_dag(range); // TODO Could return improved too

	check_consistency();

	return improved;
}

const interval var::compute_bounds() const {

	bool dummy = false;

	return lp_tighten_col(dummy);
}

void dbg_consistency(const var x[], int n) {

	assert (n > 0);
	for (int i=0; i<n; ++i) {
		x[i].check_consistency();
	}
}

const interval sum(const double c[], const var x[], int n) {

	assert(n>0);

	interval z = c[0]*x[0].bounds();

	for (int i=1; i<n; ++i) {
		z += c[i]*x[i].bounds();
	}

	return z;
}

const var lin_comb_n(const double c[], const var x[], int length) {

	dbg_consistency(x, length);

	var z(sum(c, x, length));

	if (z.bounds().degenerate()) {

		return z;
	}

	double coeff[length+1]; // TODO Maybe this can go to its own function when
	int    index[length+1]; // linear_constraint is ready

	for (int i=0; i<length; ++i) {
		coeff[i] = c[i];
		index[i] = x[i].index;
	}

	coeff[length] = -1.0;
	index[length] =  z.index;

	// 0.0 = sum coeff[i]*x[i] - 1.0*z
	var::lp->add_lin_con(0.0, coeff, index, length+1);

	z.tighten_bounds();

	return z;
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

	return bounds().contains(0.0);
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

	var::lp->add_mult_envelope(x.index, xL, xU, y.index, yL, yU, z.index);

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

	interval range(lb, ub);

	bool improved = intersect_in_dag(range);

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

	interval Z = z.bounds();

	do {

		double diam_prev = Z.diameter();

		var::lp->add_mult_envelope(y.index, Y.inf(), Y.sup(), z.index, Z.inf(), Z.sup(), x.index, improved);

		improved = z.tighten_bounds();

		using namespace std;

		cout << endl << "z: " << z << ", pass: " << ++counter << endl;

		Z = z.bounds();

		if (((diam_prev-Z.diameter()))<1.0e-6*diam_prev)
			break;

		propagate_mult(X, Y, Z);

	} while (improved);

	return z;
}

std::ostream& operator<<(std::ostream& os, const var& v) {

	return os << v.bounds() << std::flush;
}

// FIXME Eliminate this function and the example using it
void var::copy_bounds(double& lo, double& up) const {

	check_consistency();

	interval range = bounds();
	lo = range.inf();
	up = range.sup();
}

void copy_bounds(const var arr[], interval bounds[], int size) {

	for (int i=0; i<size; ++i) {

		bounds[i] = arr[i].bounds();
	}
}

void init_variables(var x[], const interval var_bounds[], int size) {

	var::reset();

	for (int i=0; i<size; ++i) {

		x[i] = var(var_bounds[i]);
	}
}

double var::width() const {

	check_consistency();

	return bounds().diameter();
}

int find_max_width(const var x[], int size) {

	assert(size>0);

	double max_width = x[0].width();

	int index = 0;

	for (int i=1; i<size; ++i) {

		const double width = x[i].width();

		if (width > max_width) {
			max_width = width;
			index = i;
		}
	}

	return index;
}

void var::reset() {

	lp->reset();
	ia_dag->reset();
}

void var::release_all() {

	delete lp;
	lp = 0;
	lp_pair::free_environment();
	delete ia_dag;
	ia_dag = 0;
}

void var::check_consistency() const {
	const interval range = bounds();
	assert(range.inf() <= range.sup());
	assert(index >= 1);
	assert(var::lp->col_type_db_or_fx(index));
}

void var::tighten_up_to(int size) {

	lp_pruning lp(var::lp, size);

	lp.prune_all();

	interval bnds[size];

	lp.copy_bounds(bnds);

	try {
		for (int i=0; i<size; ++i)
			var::ia_dag->intersect(i+1, bnds[i]);
	}
	catch (infeasible_problem& ) {
		std::cout << "Warning: numerical problems " << __FILE__ << " ";
		std::cout << __LINE__ << std::endl;
		throw numerical_problems();
	}
}

void var::tighten_all() {

	var::tighten_up_to(var::lp->n_cols());
}

}
