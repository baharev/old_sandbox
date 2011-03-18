//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010, 2011 Ali Baharev
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

//#define ASOL_DISABLE_ASSERTS

#include <ostream>
#include <algorithm>
#include "interval.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"
#include "floating_point_tol.hpp"

namespace {

void swap_if_necessary(double& lb, double& ub) {

	if (lb > ub) {
		double temp = lb;
		lb = ub;
		ub = temp;
	}
}

}

namespace asol {

interval::interval() : lb(100), ub(-100) { }

interval::interval(double value) : lb(value), ub(value) { }

interval::interval(double lo, double up) : lb(lo), ub(up) {

	ASSERT2(lb <= ub, *this);
}

interval& interval::operator+=(const interval& x) {

	ASSERT2(lb <= ub, *this);
	ASSERT2(x.lb <= x.ub, "x: "<<x);

	lb += x.lb;
	ub += x.ub;
	return *this;
}

const interval operator+(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	return interval(x.lb+y.lb, x.ub+y.ub);
}

const interval operator+(const interval& x, double y) {

	ASSERT2(x.lb <= x.ub, "x: "<<x);

	return interval(x.lb+y, x.ub+y);
}

const interval operator-(const interval& x) {

	ASSERT2(x.lb <= x.ub, "x: "<<x);

	return interval(-(x.ub), -(x.lb));
}

const interval operator-(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	return interval(x.lb-y.ub, x.ub-y.lb);
}

const interval operator-(double x, const interval& y) {

	ASSERT2(y.lb <= y.ub, "y: "<<y);

	return interval(x-y.ub, x-y.lb);
}

const interval operator*(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	double z[] = { x.lb*y.lb, x.lb*y.ub, x.ub*y.lb, x.ub*y.ub };

	double zL = *std::min_element(z, z+4);

	double zU = *std::max_element(z, z+4);

	return interval(zL, zU);
}

const interval operator*(double x, const interval& y) {

	ASSERT2(y.lb <= y.ub, "y: "<<y);

	double lb(x*y.lb), ub(x*y.ub);

	swap_if_necessary(lb, ub);

	return interval(lb, ub);
}

const interval operator*(const interval& x, double y) {

	return y*x;
}

const interval operator/(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	ASSERT2(!y.contains(0), "y: "<<y);

	double z[] = { x.lb/y.lb, x.lb/y.ub, x.ub/y.lb, x.ub/y.ub };

	double zL = *std::min_element(z, z+4);

	double zU = *std::max_element(z, z+4);

	return interval(zL, zU);
}

// Returns true and sets gap if a gap is generated, otherwise gap is undefined
bool extended_division(const interval& x, const interval& y, interval& z, interval& gap) {

	gap = interval();

	if (!y.contains(0)) {  // trivial case

		z.intersect(x/y);

		return false;
	}
	// y.contains(0)==true

	if (x.contains(0)) {  // no progress case

		return false;
	}
	// (!x.contains(0) && y.contains(0)) == true

	if (y.inf()==0 && y.sup()==0) {  // undefined case
		// FIXME Is it safe to declare it infeasible? Or should we just signal no progress?
		//ASSERT2(false, "undefined result; x, z: "<<x<<", "<<z);
		throw infeasible_problem();
	}

	return true_extended_division(x, y, z, gap);
}

bool save_gap_if_any(const double l, const double u, interval& z, interval& gap);

bool true_extended_division(const interval& x, const interval& y, interval& z, interval& gap) {

	ASSERT(!x.contains(0) && y.contains(0));

	bool ret_val = false;

	if (x.ub < 0) {

		if (y.ub==0) {

			z.prechecked_intersection(x.ub/y.lb, z.ub);
		}
		else if (y.lb==0) {

			z.prechecked_intersection(z.lb, x.ub/y.ub);
		}
		else {

			ret_val = save_gap_if_any(x.ub/y.ub, x.ub/y.lb, z, gap);
		}
	}
	else {

		if (y.ub==0) {

			z.prechecked_intersection(z.lb, x.lb/y.lb);
		}
		else if (y.lb==0) {

			z.prechecked_intersection(x.lb/y.ub, z.ub);
		}
		else {

			ret_val = save_gap_if_any(x.lb/y.lb, x.lb/y.ub, z, gap);
		}
	}

	return ret_val;
}

// FIXME It actually performs intersection which is inconsistent with intersect
bool save_gap_if_any(const double l, const double u, interval& z, interval& gap) {

	ASSERT2( l <= u, "l, u: " << l << ", " << u );

	const double zL = z.inf(), zU = z.sup();

	bool ret_val = false;

	if (zL < l && u < zU) {

		gap = interval(l, u);

		ret_val = true;
	}
	else {

		z.intersect(l, u);
	}

	return ret_val;
}

const interval sqr(const interval& x) {

	ASSERT2(x.lb <= x.ub, "x: "<<x);

	double lb(std::pow(x.lb, 2)), ub(std::pow(x.ub, 2));

	swap_if_necessary(lb, ub);

	return (x.lb<=0 && 0<=x.ub) ? interval(0, ub) : interval(lb, ub);
}

const interval sqrt(const interval& x) {

	ASSERT2(x.lb <= x.ub, "x: "<<x);
	ASSERT2(0<=x.lb, "x.lb = "<<x.lb);

	return interval(std::sqrt(x.lb), std::sqrt(x.ub));
}

const interval exp(const interval& x) {

	ASSERT2(x.lb <= x.ub, "x: "<<x);

	return interval(std::exp(x.lb), std::exp(x.ub));
}

const interval log(const interval& x) {

	ASSERT2(x.lb <= x.ub, "x: "<<x);
	ASSERT2(0<x.lb, "x.lb = "<<x.lb);

	return interval(std::log(x.lb), std::log(x.ub));
}

const interval intersection(const interval& x, const interval& y) {

	interval res(x);

	res.intersect(y);

	return res;
}

const interval hull_of(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	return interval(std::min(x.lb, y.lb), std::max(x.ub, y.ub));
}

bool interval::true_subset_of(const interval& x) const {

	ASSERT2(lb <= ub, *this);
	ASSERT2(x.lb <= x.ub, x);

	return lb >= x.lb && ub <= x.ub && (lb!=x.lb || ub !=x.ub);
}

bool interval::subset_of(const interval& x) const {

	ASSERT2(lb <= ub, *this);
	ASSERT2(x.lb <= x.ub, x);

	return lb >= x.lb && ub <= x.ub;
}

bool disjoint(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	return (x.ub < y.lb) || (y.ub < x.lb);
}

bool lessByLb(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	return x.lb < y.lb;
}

bool interval::degenerate() const {

	ASSERT2(lb <= ub, *this);

	return lb==ub;
}

bool interval::contains(double value) const {

	ASSERT2(lb <= ub, *this);

	return lb<=value && value<=ub;
}

void interval::equals(double value) {

	intersect(value, value);
}

void interval::less_than_or_equal_to(interval& rhs) {

	ASSERT2(rhs.lb <= rhs.ub, rhs);
	ASSERT2(lb <= ub, *this);

	//  this  <=  rhs
	// [a, b] <= [c, d]

	const double a = lb;

	const double d = rhs.ub;

	if (a > d) {

		throw infeasible_problem();
	}

	intersect(a, d); // b <= d; b is modified appropriately

	rhs.intersect(a, d); // a <= c; c is modified appropriately
}

bool interval::intersect(const double l, const double u) {

	ASSERT2(l <= u, "l: "<<l<<", u: "<<u);
	ASSERT2(lb <= ub, *this);

	if (is_narrow()) {
		// TODO Maybe the intersection could be computed but not written back? (May detect infeas?)
		return false;
	}

	bool improved = false;

	if (l > add_tol(lb, IMPROVEMENT_TOL)) {
		lb = l;
		improved = true;
	}

	if (u < sub_tol(ub, IMPROVEMENT_TOL)) {
		ub = u;
		improved = true;
	}

	if (lb > ub) {
		throw infeasible_problem();
	}

	return improved;
}

bool interval::prechecked_intersection(const double l, const double u) {

	ASSERT2(lb <= ub, *this);

	if (l > u) {

		throw infeasible_problem();
	}

	return intersect(l, u);
}

// z = x*y
void propagate_mult(interval& z, interval& x, interval& y) {
// TODO Eliminate propagate_mult, only used by envelopes
//	if (!x.contains(0)) { // y = z/x
//
//		y.intersect(z/x);
//	}

	interval gap;

	extended_division(z, x, y, gap);

//	if (!y.contains(0)) {  // x = z/y
//
//		x.intersect(z/y);
//	}

	extended_division(z, y, x, gap);

	// z = x*y
	z.intersect(x*y);
}

void addition_inverse(interval& z, interval& x, interval& y) {

	x.intersect(z-y);

	y.intersect(z-x);

	z.intersect(x+y);
}

void substraction_inverse(interval& z, interval& x, interval& y) {

	// z = x - y --> x = z + y
	addition_inverse(x, z, y);
}

bool division_inverse(interval& z, interval& x, interval& y, interval& gap) {

	// z = x/y --> x = z*y
	x.intersect(z*y);

	// y = x/z
	bool has_gap = extended_division(x, z, y, gap);

	z.intersect(x/y);

	return has_gap;
}

// FIXME It actually performs intersection which is inconsistent with intersect
bool sqr_inverse(interval& z, interval& x, interval& gap) {

	bool has_gap = false;

	const interval x_1 = sqrt(z);

	const interval x_2 = -x_1;

	interval x_image;

	if (disjoint(x, x_2)) {

		x_image = x_1;
	}
	else if (disjoint(x, x_1)) {

		x_image = x_2;
	}
	else {

		gap = interval(x_2.sup(), x_1.inf());

		has_gap = gap.diameter() > IMPROVEMENT_TOL; // FIXME GAP_SQRT_TOL ?

		x_image = hull_of(x_1, x_2);
	}

	x.intersect(x_image);

	z.intersect(sqr(x));

	return has_gap;
}

void exp_inverse(interval& z, interval& x) {

	x.intersect(log(z));

	z.intersect(exp(x));
}

void log_inverse(interval& z, interval& x) {

	x.intersect(exp(z));

	z.intersect(log(x));
}

// TODO Is it the best we can do?
void equality_constraint_inverse(interval& z, double rhs) {

	z.equals(rhs);
}

void copy_array(const interval src[], interval dstn[], int size) {

	for (int i=0; i<size; ++i) {

		dstn[i] = src[i];
	}
}

double interval::midpoint() const {

	ASSERT2(lb <= ub, *this);

	return (lb+ub)/2.0;
}

double interval::diameter() const {

	ASSERT2(lb <= ub, *this);

	return ub-lb;
}

double interval::radius() const {

	ASSERT2(lb <= ub, *this);

	return (ub-lb)/2.0;
}

double interval::inf() const {

	ASSERT2(lb <= ub, *this);

	return lb;
}

double interval::sup() const {

	ASSERT2(lb <= ub, *this);

	return ub;
}

bool interval::is_narrow(const double TOLERANCE) const { // FIXME Only for testing

	ASSERT2(lb <= ub, *this);

	bool ret_val = true;

	const double diameter = ub-lb;

	if (diameter >= TOLERANCE) { // diameter < TOLERANCE means narrow

		const double fabs_lb = std::fabs(lb);

		const double fabs_ub = std::fabs(ub);

		const double abs_max = fabs_lb < fabs_ub ? fabs_ub : fabs_lb;

		ASSERT(abs_max > 0);

		ret_val = (diameter/abs_max) < TOLERANCE;
	}

	return ret_val;
}

bool easy_containment(double x, const interval& y) {

	ASSERT2(y.lb <= y.ub, y);

	return sub_tol(y.lb, EASY_CONT_TOL)<=x && x<=add_tol(y.ub, EASY_CONT_TOL);
}

std::ostream& operator<<(std::ostream& os, const interval& x) {

	return os << "[ " << x.lb << ", " << x.ub << "]";
}

}
