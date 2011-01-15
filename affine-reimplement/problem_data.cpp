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
#include <functional>
#include <iterator>
#include <ostream>
#include <set>
#include "problem_data.hpp"
#include "delete_struct.hpp"
#include "diagnostics.hpp"
#include "printer.hpp"
#include "index_set.hpp"

using namespace std;

namespace {
// TODO Where should this go, same as in convert.cpp
typedef set<int> Set;

}

namespace asol {

problem_data::problem_data() {

	number_of_vars = unused_index = 0;
}

int problem_data::next_index() {

	return unused_index++;
}

int problem_data::peek_index() const {

	return unused_index;
}

void problem_data::add_variable(double lb, double ub) {

	++number_of_vars;
	initial_box.push_back(Bounds(lb, ub));;
}

void problem_data::add_primitive(Primitive* p) {

	primitives.push_back(p);
}

void problem_data::add_numeric_constant(int index, double value) {

	typedef pair<map<int,double>::iterator,bool> Result;

	Result r = numeric_constants.insert(Pair(index, value));

	ASSERT2(r.second, "index "<<index<<" already inserted");
}

void problem_data::add_common_subexpression(int index) {

	common_subexpressions.push_back(index);
}

int problem_data::last_constraint_offset() const {

	const int n = static_cast<int> ( constraints_rhs.size() );

	return n - 1;
}

int problem_data::add_constraint_rhs(double value) {

	const int primitive_index = primitives_size();

	constraints_rhs.push_back(Pair(primitive_index, value));

	return last_constraint_offset();
}

int problem_data::number_of_variables() const {

	ASSERT(number_of_vars == static_cast<int> (initial_box.size()));
	return number_of_vars;
}

const vector<primitive<builder>*>& problem_data::get_primitives() const {

	return primitives;
}

const map<int,double>& problem_data::get_numeric_constants() const {

	return numeric_constants;
}

const IntVector& problem_data::get_common_subexpressions() const {

	return common_subexpressions;
}

const PairVector& problem_data::get_rhs_of_constraints() const {

	return constraints_rhs;
}

const BoundVector& problem_data::get_initial_box() const {

	return initial_box;
}

int problem_data::primitives_size() const {

	return static_cast<int>(primitives.size());
}

problem_data::~problem_data() {

	for_each(primitives.begin(), primitives.end(), Delete());
}

void problem_data::print_info(ostream& out) const {

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


void problem_data::print_primitives(ostream& out) const {

	out << "Primitives in plain text format" << endl << endl;

	recorder* const rec = new printer(out, numeric_constants);

	const int n = primitives_size();

	for (int i=0; i<n; ++i) {

		out << i << ": ";

		primitives.at(i)->record(rec);
	}

	out << flush;

	delete rec;
}

void problem_data::print_type2_common_subexpressions(ostream& out) const {

	index_set* const rec = record_index_set();

	rec->collect_type2_common_subexpressions();

	out << "Type 2 common subexpressions:\n" << flush;

	const Set& type2_cse = rec->type2_common_subexpressions();

	copy(type2_cse.begin(), type2_cse.end(), ostream_iterator<int> (out, "\n"));

	out.flush();

	delete rec;
}

void problem_data::print_type3_common_subexpressions(ostream& out) const {

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
index_set* problem_data::record_index_set() const {

	index_set* const rec = new index_set(number_of_variables(), numeric_constants);

	for_each(primitives.begin(), primitives.end(), bind2nd(mem_fun(&Primitive::record), rec));

	rec->finished();

	return rec;
}

void problem_data::print_index_set(ostream& out) const {

	index_set* const rec = record_index_set();

	rec->print(out);

	delete rec;
}

void problem_data::common_subexpressions_type1(const int i, ostream& out) const {

	const int n = primitives_size();

	const Primitive* const p1 = primitives.at(i);

	for (int j=i+1; j<n; ++j) {

		const Primitive* const p2 = primitives.at(j);

		if (p1->common_subexpressions(p2)) {

			out << "Found in primitives: " << i << ", " << j << endl;
		}
	}
}

void problem_data::print_type1_common_subexpressions(ostream& out) const {

	out << "Type 1 common subexpressions:" << endl;

	const int n = primitives_size();

	for (int i=0; i<n; ++i) {

		common_subexpressions_type1(i, out);
	}
}

}
