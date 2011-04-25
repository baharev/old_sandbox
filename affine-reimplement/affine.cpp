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
#include "diagnostics.hpp"

namespace asol {

int affine::max_used_index(0);

std::vector<interval>* affine::v(0);

void affine::set_vector(std::vector<interval>* vec) {

	v = vec;
}

void affine::reset_counter() {

	max_used_index = 0;
}

void affine::add_noise_var(int index, double coeff) {

	noise_vars.push_back(epsilon(index, coeff));
}

void affine::make_variable() {

	ASSERT(noise_vars.empty());

	const interval rng = range();

	ASSERT2(!rng.degenerate() , "index, range: " << range_index << ", " << rng);

	add_noise_var(               0, rng.midpoint());

	add_noise_var(++max_used_index, rng.radius()  );
}

void affine::recompute_variable(int i) {

	ASSERT(noise_vars.size()==2);

	ASSERT(range_index==i);

	ASSERT(max_used_index==i);

	const interval rng = range();

	ASSERT2(!rng.degenerate() , "index, range: " << range_index << ", " << rng);

	reset_var(++max_used_index, rng);
}

void affine::reset_var(int index, const interval& rng) {

	epsilon& e_0 = noise_vars.at(0);

	ASSERT(e_0.index==0);

	e_0.coeff = rng.midpoint();

	epsilon& e_i = noise_vars.at(1);

	ASSERT(e_i.index==index);

	e_i.coeff = rng.diameter();
}

void affine::make_numeric_constant() {

	ASSERT(noise_vars.empty());

	const interval rng = range();

	ASSERT2(rng.degenerate() , "index, range: " << range_index << ", " << rng);

	add_noise_var(0, rng.inf());
}

void affine::check_if_numeric_constant() {

	ASSERT(noise_vars.size()==1);

	ASSERT2(range().degenerate() , "index, range: " << range_index << ", " << range());
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
