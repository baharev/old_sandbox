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
#include "search_procedure.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
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

	builder* box = prob->initial_box();

	prob->evaluate(box);

	delete[] box;

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

	const interval operator()(const std::pair<double,double>& p) {

		return interval(p.first, p.second);
	}
};

void search_procedure::push_initial_box_to_deque() {

	ASSERT(pending_boxes.empty());

	const BoundVector& initial_box = representation->get_initial_box();

	interval* box = new interval[n_vars];

	std::transform(initial_box.begin(), initial_box.end(), box, pair2interval());

	pending_boxes.push_back(box);
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

	ASSERT(box_orig==0);

	box_orig = pending_boxes.front();

	pending_boxes.pop_front();
}

void search_procedure::print_statistics() const {

	cout << endl;
	cout << "=========================================================" << endl;
	cout << "Number of splits: " << splits << ", solutions: ";
	cout << solutions_found << endl;
}

void search_procedure::process_box() {

	do {

		iteration_step();
	}
	while (not_discarded() && sufficient_progress());
}

void search_procedure::split_if_not_discarded() {

	if (not_discarded()) {

		split();
	}
}

bool search_procedure::not_discarded() const {

	return box_orig != 0;
}

void search_procedure::iteration_step() {

}

bool search_procedure::sufficient_progress() {

	return false;
}

void search_procedure::split() {

}

}
