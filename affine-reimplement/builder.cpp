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

#include <iostream>
#include <set>
#include <typeinfo>
#include <assert.h>
#include "builder.hpp"
#include "primitives.hpp"
#include "demangle.hpp"

namespace asol {

typedef std::set<int> Set;

int builder::unused_index = 0;

PrimVector builder::primitives = PrimVector();

PairVector builder::numeric_constants = PairVector();

IntVector builder::common_subexpressions = IntVector();

PairVector builder::constraints_rhs = PairVector();

int builder::number_of_variables() {

	return unused_index;
}

const PrimVector& builder::get_primitives() {

	return primitives;
}

const PairVector& builder::get_numeric_constants() {

	return numeric_constants;
}

const IntVector& builder::get_common_subexpressions() {

	return common_subexpressions;
}

const PairVector& builder::get_rhs_of_constraints() {

	return constraints_rhs;
}

void builder::reset() {

	unused_index = 0;

	primitives.clear();

	numeric_constants.clear();

	common_subexpressions.clear();

	constraints_rhs.clear();
}

builder::builder() : index(-1) {

}

builder::builder(double value) : index(unused_index++) {

}

builder::builder(double lb, double ub) : index(unused_index++) {

}

// FIXME Remove code duplication!
const builder operator+(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::primitives.push_back(new addition(z.index, x.index, y.index));

	return z;
}

const builder operator-(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::primitives.push_back(new substraction(z.index, x.index, y.index));

	return z;
}

const builder operator*(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::primitives.push_back(new multiplication(z.index, x.index, y.index));

	return z;
}

const builder operator/(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::primitives.push_back(new division(z.index, x.index, y.index));

	return z;
}

const builder sqr(const builder& x) {

	x.dbg_consistency();

	const builder z(0);

	builder::primitives.push_back(new square(z.index, x.index));

	return z;
}

const builder operator+(const builder& x, double y) {

	const builder Y = builder(y);

	builder::numeric_constants.push_back(Pair(Y.index, y));

	return x+Y;
}

const builder operator-(double x, const builder& y) {

	const builder X = builder(x);

	builder::numeric_constants.push_back(Pair(X.index, x));

	return X-y;
}

const builder operator*(double x, const builder& y) {

	const builder X = builder(x);

	builder::numeric_constants.push_back(Pair(X.index, x));

	return X*y;
}

// TODO Is there any benefit in saving this information?
void builder::mark_as_common_subexpression() const {

	dbg_consistency();

	common_subexpressions.push_back(index);
}

int builder::last_constraint_offset() {

	const int n = static_cast<int> ( constraints_rhs.size() );

	return n - 1;
}

void builder::equals(double value) const {

	dbg_consistency();

	const int primitive_index = static_cast<int>(primitives.size());

	constraints_rhs.push_back(Pair(primitive_index, value));

	primitives.push_back(new equality_constraint(index, last_constraint_offset()));
}

void builder::record_occurence_info() {

	const int n = static_cast<int>(constraints_rhs.size());

	for (int i=0; i<n; ++i) {

		occurence_info_of_constraint(i);
	}
}

void dump_index_set(int k, const Set& index_set) {

	std::cout << "Constraint " << k << std::endl;

	Set::const_iterator i = index_set.begin();

	while (i!=index_set.end()) {

		std::cout << *i << std::endl;

		++i;
	}
}

void builder::occurence_info_of_constraint(const int k) {

	const int start = (k!=0) ? constraints_rhs.at(k-1).first + 1 : 0;

	const int end   =          constraints_rhs.at(k  ).first;

	assert(start < end);

	Set index_set;

	for (int i=start; i<end; ++i) {

		primitives.at(i)->record_indices(index_set);
	}

	dump_index_set(k, index_set);
}

void builder::dbg_dump_type_of_primitives() {

	using namespace std;

	const int n = static_cast<int> (primitives.size());

	cout << "Dumping type of " << n << " constraints" << endl;

	for (int i=0; i<n; ++i) {

		cout << i << ": ";
		cout << demangle( typeid(*primitives.at(i)).name() ) << endl;
	}
}

void builder::dbg_consistency() const {

	assert(0<=index && index<unused_index);
}

void dbg_consistency(const builder& x, const builder& y) {

	x.dbg_consistency();
	y.dbg_consistency();
}

}
