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
#include "dag.hpp"
#include "interval.hpp"

namespace asol {

dag::dag() {

	init();
}

void dag::init() {

	interval dummy;

	variables.push_back(dummy);
}

void dag::add(int index, const interval& bounds) {

	assert(index == static_cast<int> (variables.size()));

	variables.push_back(bounds);
}

void dag::add(int index, double lb, double ub) {

	add(index, interval(lb, ub));
}

const interval dag::bounds(int index) const {

	return variables.at(index);
}

bool dag::intersect(int index, const interval& other) {

	return variables.at(index).intersect(other);
}

void dag::reset() {

	variables.clear();

	init();
}

}
