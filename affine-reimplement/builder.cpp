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

#include "builder.hpp"
#include "primitives.hpp"

namespace asol {

int builder::unused_index = 0;

std::vector<primitive*> builder::primitives = std::vector<primitive*> ();

int builder::number_of_variables() {

	return unused_index;
}

const std::vector<primitive*>& builder::get_primitives() {

	return primitives;
}

void builder::reset() {

	unused_index = 0;

	primitives.clear();
}

builder::builder() : index(-1) {

}

builder::builder(double value) : index(unused_index++) {

}

builder::builder(double lb, double ub) : index(unused_index++) {

}

const builder operator+(const builder& x, const builder& y) {

	builder z(0);

	builder::primitives.push_back(new addition(z.index, x.index, y.index));

	return z;
}

}
