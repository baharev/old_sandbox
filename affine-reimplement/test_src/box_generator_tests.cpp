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
#include <cmath>
#include <iostream>
#include <set>
#include "box_generator_tests.hpp"
#include "box_generator.hpp"
#include "diagnostics.hpp"

using namespace std;

namespace asol {

vector<interval> variables;

vector<interval> orig_box;

vector<int> index_set;

struct orderByLb : public binary_function<interval,interval,bool> {

	bool operator()(const interval& x, const interval& y) const {

		return lessByLb(x,y);
	}
};

typedef set<interval,orderByLb> Set;

vector<Set> subboxes;

void clear() {

	variables.clear();

	orig_box.clear();

	index_set.clear();

	subboxes.clear();
}

void init(const int number_of_variables) {

	clear();

	for (int i=0; i<number_of_variables; ++i) {

		const interval range = interval(i/10.0, (i+10)/10.0);

		variables.push_back(range);

		orig_box.push_back(range);

		index_set.push_back(i);
	}

	subboxes.resize(number_of_variables);
}

void register_box() {

	const int n = static_cast<int>(variables.size());

	for (int i=0; i<n; ++i) {

		subboxes.at(i).insert(variables.at(i));
	}
}

void check_counter(const int counter, const int equal_parts) {

	const int num_of_var = static_cast<int>(variables.size());

	const int expected = (int) (std::pow((double)equal_parts, num_of_var)+0.001);

	ASSERT2(counter==expected,"counter, expected: "<<counter<<", "<<expected);

	cout << "number of variables: " << num_of_var << ", " << "equal parts: ";
	cout << equal_parts << ", subboxes: " << expected << endl;
}

void check_component(const interval& component, const Set& parts, const int equal_parts) {

	Set::const_iterator part = parts.begin();

	ASSERT( part->inf() == component.inf() );

	double prev_sup = part->sup();

	++part;

	for (int i=1; i<equal_parts; ++i, ++part) {

		ASSERT( part != parts.end() );

		const double lb(part->inf()), ub(part->sup());

		ASSERT( lb < ub );

		ASSERT( prev_sup == lb );

		prev_sup = ub;
	}

	ASSERT( part == parts.end() );

	ASSERT( parts.rbegin()->sup() == component.sup() );
}

void check_coverage(const int equal_parts) {

	const int num_of_var = static_cast<int>(variables.size());

	for (int i=0; i<num_of_var; ++i) {

		const interval component = orig_box.at(i);

		const Set& parts = subboxes.at(i);

		check_component(component, parts, equal_parts);
	}
}

void run(const int equal_parts) {

	box_generator generator(variables, index_set, equal_parts);

	int loop_counter = 0;

	while (generator.get_next()) {

		variables.assign(orig_box.begin(), orig_box.end());

		generator.set_box();

		register_box();

		++loop_counter;
	}

	check_counter(loop_counter, equal_parts);

	check_coverage(equal_parts);
}

void test(const int number_of_variables, const int equal_parts) {

	init(number_of_variables);

	run(equal_parts);
}

void fail_test(const int number_of_variables, const int equal_parts) {

	try {

		test(number_of_variables, equal_parts);
	}
	catch (logic_error& e) {

		cout << "testing invalid arguments: " << e.what() << endl;
		return;
	}

	throw logic_error("testing failure failed");
}

void run_box_generator_test() {

	cout << "###############################################" << endl;
	cout << "Box generator tests" << endl;

	fail_test(0, 2);
	fail_test(2, 1);
	fail_test(2, 0);

	test(1, 2);
	test(1, 3);
	test(1, 4);
	test(2, 2);
	test(2, 3);
	test(2, 4);
	test(3, 2);
	test(3, 3);
	test(3, 4);
	test(4, 2);
	test(4, 3);
	test(4, 4);

	clear();
}

}
