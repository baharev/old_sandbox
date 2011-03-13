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

#include <functional>
#include "builder.hpp"
#include "problem_data.hpp"
#include "diagnostics.hpp"

namespace asol {

typedef addition<builder> Addition;
typedef substraction<builder> Substraction;
typedef multiplication<builder> Multiplication;
typedef division<builder> Division;
typedef square<builder> Square;
typedef exponential<builder> Exponential;
typedef logarithm<builder> Logarithm;
typedef equality_constraint<builder> Equality_constraint;
typedef common_subexpression<builder> Common_subexpression;
typedef less_than_or_equal_to<builder> Less_than_or_equal_to;

problem_data* builder::problem = new problem_data;

void builder::reset() {

	builder::release();

	problem = new problem_data;
}

void builder::release() {

	delete problem;

	problem = 0;
}

const problem_data* builder::get_problem_data() {

	return problem;
}

void builder::finished() {

	problem->build_index_set();
}

builder::builder() : index(-1) {

}

builder::builder(double ) : index(problem->next_index()) {

}

builder::builder(double lb, double ub) : index(problem->next_index()) {

	ASSERT2(lb<=ub, "lb, ub: "<<lb<<", "<<ub);

	problem->add_variable(lb, ub);
}

// FIXME This duplication is difficult to remove
const builder operator+(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::problem->add_primitive(new Addition(z.index, x.index, y.index));

	return z;
}

const builder operator-(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::problem->add_primitive(new Substraction(z.index, x.index, y.index));

	return z;
}

const builder operator*(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::problem->add_primitive(new Multiplication(z.index, x.index, y.index));

	return z;
}

const builder operator/(const builder& x, const builder& y) {

	dbg_consistency(x, y);

	const builder z(0);

	builder::problem->add_primitive(new Division(z.index, x.index, y.index));

	return z;
}

const builder sqr(const builder& x) {

	x.dbg_consistency();

	const builder z(0);

	builder::problem->add_primitive(new Square(z.index, x.index));

	return z;
}

const builder exp(const builder& x) {

	x.dbg_consistency();

	const builder z(0);

	builder::problem->add_primitive(new Exponential(z.index, x.index));

	return z;
}

const builder log(const builder& x) {

	x.dbg_consistency();

	const builder z(0);

	builder::problem->add_primitive(new Logarithm(z.index, x.index));

	return z;
}

const builder operator+(const builder& x, double y) {

	const builder Y = builder(y);

	builder::problem->add_numeric_constant(Y.index, y);

	return x+Y;
}

const builder operator+(double x, const builder& y) {

	return y+x;
}

const builder operator-(const builder& x, double y) {

	return x+(-y);
}

const builder operator-(double x, const builder& y) {

	const builder X = builder(x);

	builder::problem->add_numeric_constant(X.index, x);

	return X-y;
}

const builder operator*(double x, const builder& y) {

	const builder X = builder(x);

	builder::problem->add_numeric_constant(X.index, x);

	return X*y;
}

const builder operator/(double x, const builder& y) {

	const builder X = builder(x);

	builder::problem->add_numeric_constant(X.index, x);

	return X/y;
}

void builder::mark_as_common_subexpression() const {

	dbg_consistency();
	// FIXME It is a sort of duplication
	int ordinal = problem->add_common_subexpression(index);

	problem->add_primitive(new Common_subexpression(index, ordinal));
}

void builder::equals(double value) const {

	dbg_consistency();

	int constraint_offset = problem->add_constraint_rhs(value);

	problem->add_primitive(new Equality_constraint(index, constraint_offset, value));
}

void builder::less_than_or_equal_to(const builder& rhs) const {

	dbg_consistency();
	rhs.dbg_consistency();

	problem->add_primitive(new Less_than_or_equal_to(index, rhs.index));
}

void builder::dbg_consistency() const {

	ASSERT2(0<=index && index<problem->peek_index(), "index, unused_index: "<<index<<", "<<problem->peek_index());
}

void dbg_consistency(const builder& x, const builder& y) {

	x.dbg_consistency();
	y.dbg_consistency();
}

}
