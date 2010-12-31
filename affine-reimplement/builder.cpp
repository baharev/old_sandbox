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

#include <algorithm>
#include <iterator>
#include <ostream>
#include <set>
#include "builder.hpp"
#include "diagnostics.hpp"
#include "primitives.hpp"
#include "printer.hpp"
#include "index_set.hpp"

using namespace std;

namespace asol {

typedef set<int> Set;

int builder::unused_index = 0;

int builder::number_of_vars = 0;

PrimVector builder::primitives = PrimVector();

map<int,double> builder::numeric_constants = map<int,double> ();

IntVector builder::common_subexpressions = IntVector();

PairVector builder::constraints_rhs = PairVector();

BoundVector builder::initial_box = BoundVector();

int builder::number_of_arguments() {

	return unused_index;
}

int builder::number_of_variables() {

	ASSERT(number_of_vars == static_cast<int> (initial_box.size()))
	return number_of_vars;
}

const PrimVector& builder::get_primitives() {

	return primitives;
}

const map<int,double>& builder::get_numeric_constants() {

	return numeric_constants;
}

const IntVector& builder::get_common_subexpressions() {

	return common_subexpressions;
}

const PairVector& builder::get_rhs_of_constraints() {

	return constraints_rhs;
}

const BoundVector& builder::get_initial_box() {

	return initial_box;
}

int builder::primitives_size() {

	return static_cast<int>(primitives.size());
}

void builder::reset() {

	unused_index = 0;

	number_of_vars = 0;

	primitives.clear();

	numeric_constants.clear();

	common_subexpressions.clear();

	constraints_rhs.clear();

	initial_box.clear();
}

void builder::print_info(ostream& out) {

	ASSERT(number_of_vars == static_cast<int> (initial_box.size()));

	out << "=============================="                          << endl;
	out << "Number of variables:   " << number_of_vars               << endl;
	out << "Max used index (args): " << unused_index-1               << endl;
	out << "Common subexpressions: " << common_subexpressions.size() << endl;
	out << "Number of constraints: " << constraints_rhs.size()       << endl;
	out << "Number of primitives:  " << primitives.size()            << endl;
	out << "Numeric constants:     " << numeric_constants.size()     << endl;
	out << "=============================="                          << endl;
}

builder::builder() : index(-1) {

}

builder::builder(double ) : index(unused_index++) {

}

builder::builder(double lb, double ub) : index(unused_index++) {

	ASSERT2(lb<=ub, "lb, ub: "<<lb<<", "<<ub);

	++number_of_vars;
	initial_box.push_back(Bounds(lb, ub));
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

const builder exp(const builder& x) {

	x.dbg_consistency();

	const builder z(0);

	builder::primitives.push_back(new exponential(z.index, x.index));

	return z;
}

void builder::insert_numeric_constant(const int index, const double value) {

	typedef pair<map<int,double>::iterator,bool> Result;

	Result r = numeric_constants.insert(Pair(index, value));

	ASSERT2(r.second, "index "<<index<<" already inserted");
}

const builder operator+(const builder& x, double y) {

	const builder Y = builder(y);

	builder::insert_numeric_constant(Y.index, y);

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

	builder::insert_numeric_constant(X.index, x);

	return X-y;
}

const builder operator*(double x, const builder& y) {

	const builder X = builder(x);

	builder::insert_numeric_constant(X.index, x);

	return X*y;
}

const builder operator/(double x, const builder& y) {

	const builder X = builder(x);

	builder::insert_numeric_constant(X.index, x);

	return X/y;
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

	const int primitive_index = primitives_size();

	constraints_rhs.push_back(Pair(primitive_index, value));

	primitives.push_back(new equality_constraint(index, last_constraint_offset()));
}

void builder::print_primitives(ostream& out) {

	out << "Primitives in plain text format" << endl << endl;

	recorder* const rec = new printer(out, numeric_constants);

	const int n = primitives_size();

	for (int i=0; i<n; ++i) {

		out << i << ": ";

		primitives.at(i)->record(rec);
	}

	delete rec;
}

void builder::print_type2_common_subexpressions(ostream& out) {

	index_set* const rec = record_index_set();

	rec->collect_type2_common_subexpressions();

	out << "Type 2 common subexpressions:\n" << flush;

	const Set& type2_cse = rec->type2_common_subexpressions();

	copy(type2_cse.begin(), type2_cse.end(), ostream_iterator<int> (out, "\n"));

	out.flush();

	delete rec;
}

void builder::print_type3_common_subexpressions(ostream& out) {

	index_set* const rec = record_index_set();

	out << "Type 3 common subexpressions:\n" << flush;

	const Set& type3_cse = rec->type3_common_subexpressions();

	copy(type3_cse.begin(), type3_cse.end(), ostream_iterator<int> (out, "\n"));

	out.flush();

	delete rec;
}

// TODO Make rec member and pass it to expression_graph when ready; solves:
//      - transferring ownership
//      - eliminates duplication
index_set* builder::record_index_set() {

	index_set* const rec = new index_set(number_of_variables(), numeric_constants);

	const int n = primitives_size();

	for (int i=0; i<n; ++i) {

		primitives.at(i)->record(rec);
	}

	rec->finished();

	return rec;
}

void builder::print_index_set(ostream& out) {

	index_set* const rec = record_index_set();

	rec->print(out);

	delete rec;
}

void builder::common_subexpressions_type1(const int i, ostream& out) {

	const int n = primitives_size();

	const primitive* const p1 = primitives.at(i);

	for (int j=i+1; j<n; ++j) {

		const primitive* const p2 = primitives.at(j);

		if (p1->common_subexpressions(p2)) {

			out << "Found in primitives: " << i << ", " << j << endl;
		}
	}
}

void builder::print_type1_common_subexpressions(ostream& out) {

	out << "Type 1 common subexpressions:" << endl;

	const int n = primitives_size();

	for (int i=0; i<n; ++i) {

		common_subexpressions_type1(i, out);
	}
}

void builder::dbg_consistency() const {

	ASSERT2(0<=index && index<unused_index, "index, unused_index: "<<index<<", "<<unused_index)
}

void dbg_consistency(const builder& x, const builder& y) {

	x.dbg_consistency();
	y.dbg_consistency();
}

}
