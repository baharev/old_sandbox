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

#include <algorithm>
#include "index_recorder.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "primitives.hpp"
#include "problem_data.hpp"

using namespace std;

namespace asol {

index_recorder::index_recorder(const problem_data* prob) {

	n_vars = prob->number_of_variables();

	numeric_const = prob->get_numeric_constants();

	const vector<primitive<builder>*>& prim = prob->get_primitives();

	const int n = static_cast<int> (prim.size());

	for (int i=0; i<n; ++i) {

		pos = i;

		prim.at(i)->record(this);
	}

	// TODO Certain containers must have equal size
}

int index_recorder::primitive_index() const {

	return pos;
}

template <typename Map>
bool contains(const Map& map, const int index) {

	return map.find(index)==map.end() ? false : true;
}

bool index_recorder::is_numeric_constant(const int index) const {

	return contains(numeric_const, index);
}

bool index_recorder::is_defined_variable(const int index) const {

	return contains(def_var_indices, index);
}

bool index_recorder::is_variable(const int index) const {

	return index < n_vars;
}

void index_recorder::resolve_def_var_dependecies() {

	// TODO Start here!
	// iterate through current (start with lower n_var?)
	// if a def_var is found find in def_var_indices
	// get the dependencies from indices
	// append that to current
	// recursively resolve those def_vars too?
}

void index_recorder::push_back_current() {

	resolve_def_var_dependecies();

	indices.push_back(current);

	const int last_primitive_index = primitive_index();

	boundary.push_back(last_primitive_index);

	current.clear();

	ASSERT(indices.size() == boundary.size());
	ASSERT(constraint_end.size() <= boundary.size());
}

void index_recorder::record_unary_primitive(int z, int x) {

	ASSERT(!is_numeric_constant(z));

	record_arg(x);
}

void index_recorder::record_arg(const int index) {

	if (is_variable(index) || is_defined_variable(index)) {

		current.insert(index).second;
	}
}

void index_recorder::record_binary_primitive(int z, int x, int y) {

	record_unary_primitive(z, x);

	record_arg(y);
}

void index_recorder::addition(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_recorder::substraction(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_recorder::multiplication(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_recorder::division(int z, int x, int y) {

	record_binary_primitive(z, x, y);
}

void index_recorder::square(int z, int x) {

	record_unary_primitive(z, x);
}

void index_recorder::exponential(int z, int x) {

	record_unary_primitive(z, x);
}

void index_recorder::logarithm(int z, int x) {

	record_unary_primitive(z, x);
}

void index_recorder::equality_constraint(int , int , double ) {

	constraint_end.push_back(primitive_index());

	push_back_current();
}

void index_recorder::common_subexpression(int index, int ) {

	def_var_indices.insert(make_pair(index, indices.size()));

	push_back_current();
}

void index_recorder::less_than_or_equal_to(int lhs, int rhs) {

	ASSERT(!is_numeric_constant(rhs));

	record_unary_primitive(lhs, rhs);

	constraint_end.push_back(primitive_index());

	push_back_current();
}

}
