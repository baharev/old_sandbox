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
#include "expression_graph_test.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "expression_graph.hpp"
#include "interval.hpp"
#include "problem.hpp"
#include "problem_data.hpp"

using namespace std;

namespace asol {

const problem_data* build(const problem<builder>* prob) {

	builder::reset();

	builder* box = prob->initial_box();

	prob->evaluate(box);

	delete[] box;

	builder::finished();

	const problem_data* const representation = builder::get_problem_data();

	ASSERT(representation->number_of_variables() == prob->number_of_variables());

	delete prob;

	return representation;
}

void print_data(const problem_data* representation) {

	representation->print_info(cout);
	representation->print_primitives(cout);
	//representation->print_index_set(cout);
	representation->print_type1_common_subexpressions(cout);
	representation->print_type2_common_subexpressions(cout);
	representation->print_type3_common_subexpressions(cout);
}

void dag_test(const problem<builder>* prob) {

	const problem_data* const representation = build(prob);

	expression_graph<interval> dag(representation);

	print_data(representation);

	builder::reset();

	for (int i=0; i<40; ++i) {
		dag.revise_all();
	}

	dag.show_variables(cout);

	cout << "Last value: " << dag.last_value() << endl;
}

void print_sparsity(const problem<builder>* prob) {

	const problem_data* const representation = build(prob);

	representation->print_variable_occurences(cout);
}

void test_system_of_equations(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	for (int i=0; i<40; ++i) {
		dag.revise_all2();
	}

	dag.show_variables(cout);

	cout << "Last value: " << dag.last_value() << endl;
}

void test_iterative_revision(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	dag.iterative_revision();

	dag.show_variables(cout);

	cout << "Last value: " << dag.last_value() << endl;
}

void probing(expression_graph<interval>& dag, const interval box[], int size);

void test_probing_Jacobsen(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	interval box[] = {
			interval(0.92, 0.95),
			interval(1.0e-4, 1.0),
			interval(1.0e-4, 1.0),
			interval(1.0e-4, 1.0),
			interval(1.0e-4, 1.0),
			interval(1.0e-4, 1.0),
			interval(1.0e-4, 1.0),
			interval(1.0e-4, 1.0),
			interval(2.0, 4.0),
			interval(2.0, 4.0),
			interval(2.0, 4.0),
			interval(2.0, 4.0),
			interval(2.0, 4.0),
			interval(2.0, 4.0),
			interval(2.0, 4.0),
			interval(0.49, 0.52)
	};

	probing(dag, box, sizeof(box)/sizeof(box[0]));
}

}

