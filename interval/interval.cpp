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
#include <cmath>
#include <algorithm>
#include <assert.h>
#include "exceptions.hpp"
#include "interval.hpp"

using namespace std;

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

	assert(lb <= ub);
}

interval& interval::operator+=(const interval& x) {

	assert( (lb<=ub) && (x.lb<=x.ub) );
	lb += x.lb;
	ub += x.ub;
	return *this;
}

const interval operator+(const interval& x, const interval& y) {

	assert(x.lb <= x.ub && y.lb <= y.ub);

	return interval(x.lb+y.lb, x.ub+y.ub);
}

const interval operator+(const interval& x, double y) {

	assert(x.lb <= x.ub);

	return interval(x.lb+y, x.ub+y);
}

const interval operator-(const interval& x, const interval& y) {

	assert(x.lb <= x.ub && y.lb <= y.ub);

	return interval(x.lb-y.ub, x.ub-y.lb);
}

const interval operator-(double x, const interval& y) {

	assert(y.lb <= y.ub);

	return interval(x-y.ub, x-y.lb);
}

const interval operator*(const interval& x, const interval& y) {

	assert(x.lb <= x.ub && y.lb <= y.ub);

	double z[] = { x.lb*y.lb, x.lb*y.ub, x.ub*y.lb, x.ub*y.ub };

	double zL = *min_element(z, z+4);

	double zU = *max_element(z, z+4);

	return interval(zL, zU);
}

const interval operator*(double x, const interval& y) {

	assert(y.lb <= y.ub);

	double lb(x*y.lb), ub(x*y.ub);

	swap_if_necessary(lb, ub);

	return interval(lb, ub);
}

const interval operator/(const interval& x, const interval& y) {

	assert(x.lb <= x.ub && y.lb <= y.ub);

	if (y.contains(0)) {

		throw assertion_error();
	}

	double z[] = { x.lb/y.lb, x.lb/y.ub, x.ub/y.lb, x.ub/y.ub };

	double zL = *min_element(z, z+4);

	double zU = *max_element(z, z+4);

	return interval(zL, zU);
}

const interval sqr(const interval& x) {

	assert(x.lb <= x.ub);

	double lb(std::pow(x.lb, 2)), ub(std::pow(x.ub, 2));

	swap_if_necessary(lb, ub);

	return (x.lb<=0 && 0<=x.ub) ? interval(0, ub) : interval(lb, ub);
}

bool interval::degenerate() const {

	assert(lb <= ub);

	return lb==ub;
}

bool interval::contains(double value) const {

	assert(lb <= ub);

	return lb<=value && value<=ub;
}

bool interval::intersect(const interval& other) {

	return intersect(other.inf(), other.sup());
}

bool interval::intersect(const double l, const double u) {

	assert(l<=u);
	assert(lb<=ub);

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


void copy_array(const interval src[], interval dstn[], int size) {

	for (int i=0; i<size; ++i) {

		dstn[i] = src[i];
	}
}

double interval::midpoint() const {

	assert(lb <= ub);

	return (lb+ub)/2.0;
}

double interval::diameter() const {

	assert(lb <= ub);

	return ub-lb;
}

double interval::radius() const {

	assert(lb <= ub);

	return (ub-lb)/2.0;
}

double interval::inf() const {

	assert(lb <= ub);

	return lb;
}

double interval::sup() const {

	assert (lb <= ub);

	return ub;
}

std::ostream& operator<<(std::ostream& os, const interval& x) {

	return os << "[ " << x.inf() << ", " << x.sup() << "]";
}

}
