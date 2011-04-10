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
#include <iostream>
#include <iterator>
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

	ASSERT2(current.empty(),"problem must be closed by a constraint");

	compute_constraint_index_set();

	// TODO Check if all variables and def_vars are really contained in the index set
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

const set<int> index_recorder::extract_variables(const set<int>& s) {

	return set<int>(s.begin(), s.lower_bound(n_vars));
}

void index_recorder::resolve_def_var_dependecies() {

	typedef set<int>::const_iterator itr;

	itr i = current.lower_bound(n_vars);

	for ( ; i!= current.end(); ++i)	{

		map<int,int>::const_iterator k = def_var_indices.find(*i);

		ASSERT2(k!=def_var_indices.end(), "index not found: " << (*i) );

		const set<int>& dependencies = indices.at(k->second);

		const set<int> vars = extract_variables(dependencies); // recursively resolve def_vars ?

		current.insert(vars.begin(), vars.end());
	}
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

void index_recorder::equality_constraint(int , int , double ) {

	constraint_end.push_back(primitive_index());

	push_back_current();
}

void index_recorder::common_subexpression(int index, int ) {

	def_var_indices.insert(make_pair(index, indices.size()));

	push_back_current();
}

void index_recorder::less_than_or_equal_to(int lhs, int rhs) {

	ASSERT(!is_numeric_constant(lhs) && !is_numeric_constant(rhs));

	record_arg(lhs);

	record_arg(rhs);

	constraint_end.push_back(primitive_index());

	push_back_current();
}

void index_recorder::record_arg(const int index) {

	if (is_variable(index) || is_defined_variable(index)) {

		current.insert(index).second;
	}
}

template <typename T>
void dump(const vector<T>& index, const vector<int>& last_primitive) {

	const int n = static_cast<int>(index.size());

	for (int i=0; i<n; ++i) {

		cout << i << ": ";

		const T& s = index.at(i);

		copy(s.begin(), s.end(), ostream_iterator<int>(cout, "\t"));

		cout << '(' << last_primitive.at(i) << ')' << endl;
	}

	cout << endl;
}

void index_recorder::dump() const {
	// TODO Ask on SO why asol:: is needed
	asol::dump(indices, boundary);

	asol::dump(constraint_indices, constraint_end);
}

void index_recorder::merge_up_to(const int end) {

	ASSERT(current.empty());

	do {

		++idx;

		const set<int>& s = indices.at(idx);

		current.insert(s.begin(), s.end());
	}
	while (boundary.at(idx)!=end);

	ASSERT(!current.empty());

	// TODO Only variables are pushed back, how about defined vars?
	set<int> vars = extract_variables(current);

	constraint_indices.push_back(vector<int>(vars.begin(), vars.end()));

	current.clear();
}

void index_recorder::compute_constraint_index_set() {

	idx = -1;

	const int n = constraint_end.size();

	for (int i=0; i<n; ++i) {

		merge_up_to(constraint_end.at(i));
	}
}

const vector<vector<int> >& index_recorder::constraint_index_sets() const {

	return constraint_indices;
}

void index_recorder::record_unary_primitive(int z, int x) {

	ASSERT(!is_numeric_constant(z));

	record_arg(x);
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

}
