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

#include <ostream>
#include <cmath>
#include <algorithm>
#include "interval.hpp"
#include "exceptions.hpp"
#include "diagnostics.hpp"

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

const interval operator/(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	ASSERT2(!y.contains(0), "y: "<<y);

	double z[] = { x.lb/y.lb, x.lb/y.ub, x.ub/y.lb, x.ub/y.ub };

	double zL = *std::min_element(z, z+4);

	double zU = *std::max_element(z, z+4);

	return interval(zL, zU);
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

const interval hull_of(const interval& x, const interval& y) {

	ASSERT2(x.lb<=x.ub && y.lb<=y.ub, "x: "<<x<<", y: "<<y);

	return interval(std::min(x.lb, y.lb), std::max(x.ub, y.ub));
}

bool interval::subset_of(const interval& x) const {

	ASSERT2(lb <= ub, *this);
	ASSERT2(x.lb <= x.ub, x);

	return lb >= x.lb && ub <= x.ub && (lb!=x.lb || ub !=x.ub);
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

bool interval::intersect(const interval& other) {

	return intersect(other.inf(), other.sup());
}

void interval::assign(const interval& other) {

	intersect(other);
}

void interval::equals(double value) {

	intersect(value, value);
}

bool interval::intersect(const double l, const double u) {

	ASSERT2(l <= u, "l: "<<l<<", u: "<<u);
	ASSERT2(lb <= ub, *this);

	bool improved = false;

	if (l>lb) {
		lb = l;
		improved = true;
	}

	if (u<ub) {
		ub = u;
		improved = true;
	}

	if (lb>ub) {
		throw infeasible_problem();
	}

	return improved;
}

// z = x*y
void propagate_mult(interval& z, interval& x, interval& y) {

	if (!x.contains(0)) { // y = z/x

		y.intersect(z/x);
	}

	if (!y.contains(0)) {  // x = z/y

		x.intersect(z/y);
	}

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

void multiplication_inverse(interval& z, interval& x, interval& y) {
	// TODO Use extended division and save the gap?
	// TODO Eliminate propagate_mult, only used by envelopes
	propagate_mult(z, x, y);
}

void division_inverse(interval& z, interval& x, interval& y) {

	ASSERT2(!y.contains(0), "y: "<<y);

	// z = x/y --> x = z*y
	multiplication_inverse(x, z, y);
}

void sqr_inverse(interval& z, interval& x) {

	interval x_new = sqrt(z);

	if (x.inf()>=0) {
		;
	}
	else if (x.sup()<=0) {
		x_new = -x_new;
	}
	else {
		x_new = interval(-(x_new.sup()), x_new.sup());
	}

	x.intersect(x_new);

	z.intersect(sqr(x));
}

void exp_inverse(interval& z, interval& x) {

	x.intersect(log(z));

	z.intersect(exp(x));
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

bool is_narrow(const interval& x) {

	return x.diameter() < NARROW;
}

std::ostream& operator<<(std::ostream& os, const interval& x) {

	return os << "[ " << x.lb << ", " << x.ub << "]";
}

}
