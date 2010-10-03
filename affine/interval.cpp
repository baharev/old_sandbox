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

#include <iostream>
#include <assert.h>
#include "interval.hpp"

using namespace std;

namespace asol {

interval::interval() : lb(100), ub(-100) { }

interval::interval(double lo, double up) : lb(lo), ub(up) {

	assert (lb <= ub);
}

const interval operator+(const interval& x, const interval& y) {
	assert(x.lb <= x.ub && y.lb <= y.ub);

	return interval(x.lb + y.lb, x.ub + y.ub);
}

double interval::midpoint() const {

	assert (lb <= ub);

	return (lb+ub)/2.0;
}

double interval::radius() const {

	assert (lb <= ub);

	return (ub-lb)/2.0;
}

double interval::inf() const {

	assert (lb <= ub);

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