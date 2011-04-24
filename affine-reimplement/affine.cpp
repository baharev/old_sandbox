//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
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
#include "affine.hpp"

namespace asol {

int affine::max_used_index(0);

std::vector<interval>* affine::v(0);

affine::affine() : range_index(-1) {

}

// TODO Figure out how the affine expr_graph will be initialized
affine::affine(double value) {
	// Range index?
}

affine::affine(double lb, double ub) {

}

std::ostream& operator<<(std::ostream& os, const affine& x) {

	for (int i=0; i<x.size(); ++i) {

		epsilon e = x.noise_vars.at(i);

		os << e.index << ": " << e.coeff << '\n';
	}

	os << "IA: " << x.range() << "\n\n";

	return os;
}

const affine operator+(const affine& x, const affine& y) {

	return affine();
}

const affine operator-(const affine& x, const affine& y) {

	return affine();
}

const affine operator*(const affine& x, const affine& y) {

	return affine();
}

const affine operator/(const affine& x, const affine& y) {

	return affine();
}

void affine::assign(const affine& other) {

}

void affine::equals(double value) {

}

void affine::less_than_or_equal_to(affine& rhs) {

}

const affine exp(const affine& x) {

	return affine();
}

const affine log(const affine& x) {

	return affine();
}

const affine sqr(const affine& x) {

	return affine();
}

affine::~affine() {
	// Out of line dtor to make the compiler shut up
}

}
