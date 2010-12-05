//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010 Ali Baharev
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
#include <cmath> // FIXME
#include <assert.h>
#include "algorithm.hpp"
#include "constants.hpp"
#include "envelope.hpp"
#include "exceptions.hpp"
#include "problem.hpp"

namespace asol {

algorithm::algorithm(const problem* const p)
	: n(p->size()), prob(p), box(new var[n]), box_orig(0)
{
	solutions_found = splits = depth = boxes_processed = 0;
}

algorithm::~algorithm() {

	assert(pending.empty()); // TODO Find a better way

	assert(box_orig==0);

	delete[] box;

	var::release_all();
}

void algorithm::run() {

	add_initial_box();

	do {

		set_current_box();

		iteration_step();

	} while (!pending.empty());

	print_statistics();
}

void algorithm::add_initial_box() {

	pending.push_back(prob->initial_box());
}

void algorithm::set_current_box() {

	assert(box_orig==0);

	box_orig = pending.front();

	pending.pop_front();

	init_variables(box, box_orig, n);
}

void algorithm::prepare_to_repeat() {

	std::cout << "Repeating pruning steps" << std::endl;

	pending.push_front(box_orig); // FIXME Assumes box -> box_orig copy

	box_orig = 0;
}

void algorithm::iteration_step() {

	try {
		contracting_step();
	}
	catch (infeasible_problem& ) {
		delete_box();
		return;
	}
	catch (numerical_problems& ) {
		rollback();
	}
	catch (convergence_reached& ) {
		print_box();
		delete_box();
		return;
	}

	if (sufficient_progress()) {
		prepare_to_repeat();
	}
	else {
		split();
	}
}

void algorithm::contracting_step() {

	increment_counters();

	evaluate();

	try {

		lp_pruning();
	}
	catch (infeasible_problem& ) {

		std::cout << "Warning: numerical problems in LP pruning" << std::endl;

		throw numerical_problems();
	}
}

void algorithm::increment_counters() {

	// TODO Track and log depth info too
	++boxes_processed;
}

void algorithm::evaluate() {

	std::cout << "Evaluation" << std::endl;

	prob->evaluate(box);

	check_convergence();
}

// TODO Replace this mock implementation with Acterberg's heuristic
void algorithm::lp_pruning() {

	std::cout << "Running LP pruning" << std::endl;
	// FIXME Check for convergence before tightening a variable
	for (int i=0; i<n; ++i) {

		box[i].tighten_bounds();
	}

	check_convergence();
}

void algorithm::delete_box() {

	std::cout << "Box discarded" << std::endl;

	delete[] box_orig;

	box_orig = 0;
}

void algorithm::rollback() {

	std::cout << "Warning: numerical problems detected" << std::endl;

	init_variables(box, box_orig, n);
}

void algorithm::print_box() const {

	for (int i=0; i<n; ++i) {

		std::cout << box[i] << std::endl;
	}
}

void algorithm::check_convergence() {

	int index = find_max_width(box, n);

	if (box[index].width() <= TOL_SOLVED) {

		std::cout << "Found a solution!" << std::endl;

		++solutions_found;

		throw convergence_reached();
	}
}

int algorithm::select_index_to_split() const {

	double x1 = box[0].width();
	double D  = box[15].width();

	int index = (x1 > D)? 0 : 15;

	//int index = find_max_width(box, n);

	std::cout << "Splitting " << index << ", " << box_orig[index] << std::endl;
	std::cout << "var[index] = " << box[index] << std::endl;

	assert ( std::fabs(box[index].width()- box_orig[index].diameter()) < 1.0e-6); // FIXME
	assert ( box[index].width() > TOL_SOLVED);

	return index;
}

void algorithm::split() {

	interval* const box_new = new interval[n];

	copy_array(box_orig, box_new, n);

	const int index = select_index_to_split();

	double lb  = box_orig[index].inf();
	double ub  = box_orig[index].sup();
	double mid = box_orig[index].midpoint();

	box_orig[index] = interval(lb, mid);
	box_new[index]  = interval(mid, ub);

	pending.push_back(box_orig);
	pending.push_back(box_new);
	++splits;

	box_orig = 0;
}

// TODO Replace mock implementation
bool algorithm::sufficient_progress() {

	std::cout << "Computing progress" << std::endl;

	interval box_contracted[n];

	copy_bounds(box, box_contracted, n);

	bool sufficient = false;

	double max_reduction = 10;

	for (int i=0; i<n; ++i) {

		if (box_orig[i].diameter() == 0) {
			continue;
		}

		double reduction = box_contracted[i].diameter() / box_orig[i].diameter();

		if (reduction<0.75) {
			sufficient = true;
		}

		if (reduction<max_reduction) {
			max_reduction = reduction;
		}
	}

	std::cout << "Max progress: " << 1.0-max_reduction << std::endl;

	copy_bounds(box, box_orig, n);

	return sufficient;
}

void algorithm::print_statistics() const {

	std::cout << "Number of splits: " << splits << ", solutions: ";
	std::cout << solutions_found << std::endl;
}

}

