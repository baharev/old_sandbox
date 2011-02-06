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
#include <iostream>
#include "expression_graph_test.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "expression_graph.hpp"
#include "interval.hpp"
#include "problem.hpp"
#include "problem_data.hpp"

using namespace std;

namespace {

vector<vector<asol::interval> > sol_boxes;

vector<vector<double> > sol_vectors;

const double AMOUNT(0.01);

}

namespace asol {

void evaluate(const problem<builder>* prob) {

	builder::reset();

	builder* box = prob->initial_box();

	prob->evaluate(box);

	delete[] box;

	builder::finished();
}

class inflation {

	const interval infl;

public:

	inflation(double amount) : infl(1-amount, 1+amount) { }

	const interval operator()(const pair<double,double>& orig, const double sol) const {

		return intersection(interval(orig.first, orig.second), infl*sol);
	}
};

void copy_solution_vectors(const problem<builder>* prob) {

	const int n_var = prob->number_of_variables();

	const int n_sol = prob->number_of_stored_solutions();

	sol_vectors.resize(n_sol);

	for (int i=0; i<n_sol; ++i) {

		const double* const x = prob->solution(i);

		sol_vectors.at(i).assign(x, x + n_var);
	}
}

void copy_solutions(const problem<builder>* prob) {

	copy_solution_vectors(prob);

	const BoundVector& initial_box = builder::get_problem_data()->get_initial_box();

	const int n_sol = prob->number_of_stored_solutions();

	sol_boxes.resize(n_sol);

	for (int i=0; i<n_sol; ++i) { // Eliminating the loop would make the intentions obscure

		const double* const x = prob->solution(i);

		vector<interval>& sol= sol_boxes.at(i);

		sol.resize(initial_box.size());

		transform(initial_box.begin(), initial_box.end(), x, sol.begin(), inflation(AMOUNT));
	}
}

const problem_data* build(const problem<builder>* prob) {

	evaluate(prob);

	copy_solutions(prob);

	const problem_data* const representation = builder::get_problem_data();

	ASSERT(representation->number_of_variables() == prob->number_of_variables());

	delete prob; // FIXME Ownership is taken, however auto_ptr would make the code messy

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

	builder::reset();
}

void test_system_of_equations(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob)); // FIXME Duplication! auto_ptr?

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

void probing(expression_graph<interval>& dag);

void test_Bratu(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	probing(dag);
}

void test_Bratu_solutions(expression_graph<interval>& dag, const int n_sol) {

	for (int i=0; i<n_sol; ++i) {

		cout << endl << "Testing solution " << (i+1) << " of " << n_sol << endl;

		vector<interval>& box = sol_boxes.at(i);

		dag.set_box(&(box.at(0)), box.size());

		dag.revise_all();

		dag.show_variables(cout);

		ASSERT(dag.contains(sol_vectors.at(i)) == STRICT_CONTAINMENT);
	}
}

void test_Bratu_solutions(const problem<builder>* prob) {

	const int n_sol = prob->number_of_stored_solutions();

	expression_graph<interval> dag(build(prob));

	builder::reset();

	ASSERT(n_sol==static_cast<int>(sol_boxes.size()));

	test_Bratu_solutions(dag, n_sol);
}

}

