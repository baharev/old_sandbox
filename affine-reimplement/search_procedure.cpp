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
#include <cmath>
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
using std::fabs;

namespace asol {

search_procedure::search_procedure(const problem<builder>* p)
: prob(p), n_vars(prob->number_of_variables()), representation(0), box_orig(0)
{
	build_problem_representation();

	expr_graph = new expression_graph<interval>(representation, prob->solutions());

	push_initial_box_to_deque();

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
		// FIXME We should dump the previous one...
		if(expr_graph->contains_solution()) {
			expr_graph->dump();
			ASSERT(false);
		}

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

	expr_graph->save_containment_info();

	//expr_graph->probing();

	expr_graph->iterative_revision();

	expr_graph->check_transitions_since_last_call();

	check_convergence();
}

const double CONVERGENCE_TOL = 0.05; // FIXME Just for testing

struct wide {

	bool operator()(const interval& x) const { return !x.is_narrow(CONVERGENCE_TOL); }
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

	cout << "Box discarded" << endl; // TODO Somewhat misplaced for true solutions

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

struct diam_reduction {

	double operator()(const interval& x, const interval& y) {

		ASSERT(x.subset_of(y));

		const double y_diam = y.diameter();

		return (y_diam == 0) ? 10 : x.diameter() / y_diam; // TODO Magic number
	}
};

double search_procedure::compute_max_progress() const {

	const interval* const box_contracted = expr_graph->get_box();

	double reduction[n_vars];

	std::transform(box_contracted, box_contracted+n_vars, box_orig, reduction, diam_reduction());

	const double best_reduction = *std::min_element(reduction, reduction+n_vars);

	cout << "Best reduction: " << best_reduction << endl;

	// Check convergence was already called
	ASSERT2(best_reduction!=10,"all components have zero width"); // FIXME Magic number

	return best_reduction;
}

void search_procedure::split() {

	interval* const box_new = new interval[n_vars];

	std::copy(box_orig, box_orig+n_vars, box_new);

	//double x1 = box_orig[0].diameter();
	//double D  = box_orig[15].diameter();

	//int index = (x1 > D)? 0 : 15;

	const int index = select_index_to_split();

	double lb  = box_orig[index].inf();
	double ub  = box_orig[index].sup();

	double mid;

	if (lb <= 0 && 0 <= ub) {
		mid = (fabs(lb) > fabs(ub)) ? (lb/10.1) : (ub/10.1);
	}
	else {
		mid = box_orig[index].midpoint();
	}

	box_orig[index] = interval(lb, mid);
	box_new[index]  = interval(mid, ub);

	pending_boxes.push_back(box_orig);
	pending_boxes.push_back(box_new);

	++splits;

	box_orig = 0;
}

struct width {

	double operator()(const interval& x) { return x.diameter(); } // TODO rel diam?
};

int search_procedure::select_index_to_split() const {

	double reldiam[n_vars];

	std::transform(box_orig, box_orig+n_vars, reldiam, width());

	int index = std::max_element(reldiam, reldiam+n_vars) - reldiam;

	cout << "Splitting " << index << ", " << box_orig[index] << endl;

	//ASSERT ( ! box_orig[index].is_narrow(CONVERGENCE_TOL) ); // TODO Inconsistent with width

	return index;
}

}
