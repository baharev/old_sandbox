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
#include <iomanip>
#include "expression_graph_test.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "expression_graph.hpp"
#include "floating_point_tol.hpp"
#include "gap_probing.hpp"
#include "interval.hpp"
#include "problem.hpp"
#include "problem_data.hpp"
#include "sol_tracker.hpp"

using namespace std;

namespace {

typedef vector<asol::interval> ivector;

typedef vector<double> dvector;

vector<ivector> sol_boxes;

vector<dvector> sol_vectors;

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

	const double amount;

public:

	inflation(double amount) : amount(amount) { }

	const interval operator()(const pair<double,double>& orig, const double sol) const {

		const double l = orig.first;
		const double u = orig.second;

		ASSERT2(l<=sol&&sol<=u,"l, u, sol: "<<l<<", "<<u<<", "<<sol);

		const double lb = sub_tol(sol, amount);
		const double ub = add_tol(sol, amount);

		return intersection(interval(l, u), interval(lb, ub));
	}
};

void copy_solution_vectors(const problem<builder>* prob) {

	const int n_var = prob->number_of_variables();

	const int n_sol = prob->number_of_stored_solutions();

	sol_vectors.resize(n_sol);

	const DoubleArray2D solution_vectors = prob->solutions();

	for (int i=0; i<n_sol; ++i) {

		const double* const x = &solution_vectors[i][0];

		sol_vectors.at(i).assign(x, x + n_var);
	}
}

void copy_solutions(const problem<builder>* prob) {

	copy_solution_vectors(prob);

	const BoundVector& initial_box = builder::get_problem_data()->get_initial_box();

	const int n_sol = prob->number_of_stored_solutions();

	sol_boxes.resize(n_sol);

	const DoubleArray2D solution_vectors = prob->solutions();

	for (int i=0; i<n_sol; ++i) { // Eliminating the loop would make the intentions obscure

		const double* const x = &solution_vectors[i][0];

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
	//representation->print_index_set(cout); // FIXME Reconsider this kind of index set!
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

// FIXME Reconsider and perhaps remove this kind of index set
//void print_index_set(const problem<builder>* prob) {
//
//	const problem_data* const representation = build(prob);
//
//	representation->print_index_set(cout);
//
//	builder::reset();
//}

void print_sparsity(const problem<builder>* prob) {

	const problem_data* const representation = build(prob);

	representation->print_variable_occurences(cout);

	builder::reset();
}

template <typename MemFun>
void test_solutions(expression_graph<interval>& dag, MemFun f) {

	const int n_sol = static_cast<int> (sol_boxes.size());

	for (int i=0; i<n_sol; ++i) {

		cout << endl << "Testing solution " << (i+1) << " of " << n_sol << endl;

		ivector& box = sol_boxes.at(i);

		dag.set_box(&(box.at(0)), box.size());

		sol_tracker tracker(sol_vectors);

		tracker.save_containment_info(dag.get_box());

		(dag.*f)();

		dag.show_variables(cout);

		tracker.check_transitions_since_last_call(dag.get_box());
	}
}

void test_solutions_revise(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob)); // FIXME Duplication!

	builder::reset();

	test_solutions(dag, &expression_graph<interval>::revise_all);
}

void test_solutions_revise2(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	test_solutions(dag, &expression_graph<interval>::revise_all2);
}

void test_solutions_iterative_revise(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();
	// TODO Decide whether the gaps should be saved
	test_solutions(dag, &expression_graph<interval>::iterative_revision_save_gaps);
}

void test_solutions_probing(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	test_solutions(dag, &expression_graph<interval>::probing);
}

void test_probing_on_initial_box(const problem<builder>* prob) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	dag.probing();

	dag.show_variables(cout);
}

void extended_division_test(const problem<builder>* prob, const interval* box, const double* sol, int length) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	dag.set_box(box, length);

	dvector solution = dvector(sol, sol+length);

	cout << endl << "Initial box:" << endl;

	sol_tracker tracker(sol_vectors);

	tracker.save_containment_info(dag.get_box());

	dag.show_variables(cout);

	dag.probing();

	cout << endl << "After probing:" << endl;

	dag.show_variables(cout);

	tracker.check_transitions_since_last_call(dag.get_box());

}

void gap_probing_test(const problem<builder>* prob, interval* box, const double* sol, int length) {

	expression_graph<interval> dag(build(prob));

	builder::reset();

	dag.set_box(box, length);

	sol_tracker tracker(sol_vectors);

	tracker.save_containment_info(dag.get_box());

	dvector solution = dvector(sol, sol+length);

	cout << endl << "Initial box:" << endl;

	dag.show_variables(cout);

	interval* initial_box = new interval[length];

	copy(box, box+length, initial_box);

	gap_probing p(dag, initial_box, length);

	interval* reduced_box = p.contracted_box();

	cout << endl << "After gap probing:" << endl;

	dag.show_variables(cout);

	tracker.check_transitions_since_last_call(dag.get_box());

	delete[] reduced_box;
}

}

