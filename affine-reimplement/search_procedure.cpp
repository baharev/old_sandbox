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
#include "search_procedure.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"
#include "expression_graph.hpp"
#include "interval.hpp"
#include "problem.hpp"
#include "problem_data.hpp"

using std::cout;
using std::endl;

namespace asol {

search_procedure::search_procedure(const problem<builder>* p)
: n_vars(prob->number_of_variables()), prob(p), representation(0), box_orig(0)
{
	build_problem_representation();

	expr_graph = new expression_graph<interval>(representation);

	push_initial_box_to_deque();

	copy_stored_solutions_if_any();

	solutions_found = splits = boxes_processed = 0;

	representation = 0;

	delete prob;

	prob = 0;

	builder::reset();
}

search_procedure::~search_procedure() {

	delete expr_graph;
}

void search_procedure::build_problem_representation() {

	evaluate_with_builder();

	representation = builder::get_problem_data();

	ASSERT(representation->number_of_variables() == prob->number_of_variables());
}

void search_procedure::evaluate_with_builder() const {

	builder::reset();

	builder* x = prob->initial_box();

	prob->evaluate(x);

	delete[] x;

	builder::finished();
}

void search_procedure::copy_stored_solutions_if_any() {

	const int n_sol = prob->number_of_stored_solutions();

	sol_vectors.resize(n_sol);

	for (int i=0; i<n_sol; ++i) {

		const double* const x = prob->solution(i);

		sol_vectors.at(i).assign(x, x + n_vars);
	}
}

struct pair2interval {

	const interval operator()(const std::pair<double,double>& p) const {

		return interval(p.first, p.second);
	}
};

void search_procedure::push_initial_box_to_deque() {

	ASSERT(pending_boxes.empty());

	const BoundVector& initial_box = representation->get_initial_box();

	interval* x = new interval[n_vars];

	std::transform(initial_box.begin(), initial_box.end(), x, pair2interval());

	pending_boxes.push_back(x);
}

void search_procedure::run() {

	while (has_more_boxes()) {

		get_next_box();

		process_box();

		split_if_not_discarded();
	}

	print_statistics();
}

bool search_procedure::has_more_boxes() const {

	return !pending_boxes.empty();
}

void search_procedure::get_next_box() {

	cout << "=========================================================" << endl;

	ASSERT(box_orig == 0);

	box_orig = pending_boxes.front();

	pending_boxes.pop_front();
}

void search_procedure::print_statistics() const {

	ASSERT(2*splits+1 == boxes_processed);

	cout << endl;
	cout << "=========================================================" << endl;
	cout << "Number of splits: " << splits << ", solutions: ";
	cout << solutions_found << endl;
}

void search_procedure::process_box() {

	do {

		iteration_step();
	}
	while (not_done_with_box() && sufficient_progress());

	++boxes_processed;
}

void search_procedure::split_if_not_discarded() {

	if (not_done_with_box()) {

		split();
	}
}

bool search_procedure::not_done_with_box() const {

	return box_orig != 0;
}

void search_procedure::iteration_step() {

	try {

		contracting_step();
	}
	catch (infeasible_problem& ) {

		delete_box();
	}
	catch (numerical_problems& ) {

		ASSERT2(false,"implementation not updated properly");
	}
	catch (convergence_reached& ) {

		print_box();

		delete_box();
	}
}

void search_procedure::contracting_step() {

	expr_graph->set_box(box_orig, n_vars);

	expr_graph->probing();

	check_convergence();
}

struct wide {

	bool operator()(const interval& x) const { return !x.is_narrow(0.05); }
};

void search_procedure::check_convergence() {

	// TODO Move convergence check to expression_graph?
	const interval* const box = expr_graph->get_box();

	const interval* const elem = std::find_if(box, box+n_vars, wide());

	if (elem == box+n_vars) {

		cout << "Found a solution!" << endl;

		++solutions_found;

		throw convergence_reached();
	}
}

void search_procedure::delete_box() {

	cout << "Box discarded" << endl;

	delete[] box_orig;

	box_orig = 0;
}

void search_procedure::print_box() const {

	expr_graph->show_variables(cout);
}

bool search_procedure::sufficient(const double max_progress) const {

	return max_progress < 0.75; // TODO Introduce options class
}

bool search_procedure::sufficient_progress() {

	cout << "Computing progress" << endl;

	const double best_reduction = compute_max_progress();

	if (sufficient(best_reduction)) {

		cout << "Sufficient progress made" << endl;
		cout << "-----------------------------------------------------" << endl;
	}

	const interval* const box = expr_graph->get_box();

	std::copy(box, box+n_vars, box_orig);

	return sufficient(best_reduction);
}

struct reduction {

	double operator()(const interval& x, const interval& y) {

		// FIXME Continue from here!!!
	}
};

double search_procedure::compute_max_progress() const {

	const interval* const box_contracted = expr_graph->get_box();

	double best_reduction = 10;

	for (int i=0; i<n_vars; ++i) {

		const double orig_diam = box_orig[i].diameter();

		if (orig_diam == 0) {

			cout << "Warning: index " << i << " has zero width" << endl;

			continue;
		}

		// TODO Check if orig contains contracted

		double reduction = box_contracted[i].diameter() / orig_diam;

		if (reduction<best_reduction) {

			best_reduction = reduction;
		}
	}

	cout << "Best reduction: " << best_reduction << endl;

	return best_reduction;
}

void search_procedure::split() {

}

}
