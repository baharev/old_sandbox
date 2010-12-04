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
#include <assert.h>
#include "algorithm.hpp"
#include "envelope.hpp"
#include "exceptions.hpp"
#include "problem.hpp"

namespace asol {

algorithm::algorithm(const problem* const p)
	: n(p->size()), prob(p), box(new var[n]), box_orig(0)
{
	depth = boxes_processed = 0;
}

algorithm::~algorithm() {

	assert(pending.empty()); // TODO Find a better way

	assert(!box_orig);

	delete[] box;

	var::release_all();
}

void algorithm::run() {

	add_initial_box();

	do {

		set_current_box();

		iteration_step();

	} while (!pending.empty());
}

void algorithm::add_initial_box() {

	pending.push_back(prob->initial_box());
}

void algorithm::set_current_box() {

	box_orig = pending.front();

	pending.pop_front();

	init_variables(box, box_orig, n);
}

void algorithm::prepare_to_repeat() {

	pending.push_front(box_orig); // FIXME Assumes box -> box_orig copy
}

void algorithm::iteration_step() {

	try {

		contracting_step();
	}
	catch (infeasible_problem& ) {

		std::cout << "Box discarded" << std::endl;

		delete_box();

		return;
	}
	catch (numerical_problems& ) {

		std::cout << "Warning: numerical problems detected" << std::endl;

		rollback();
	}

	// TODO Check for convergence

	if (sufficient_progress()) {

		std::cout << "Repeating pruning steps" << std::endl;

		prepare_to_repeat();

		return;
	}

	split();
}

void algorithm::contracting_step() {

	increment_counters();

	evaluate();

	// FIXME Check for convergence here!

	std::cout << "Running LP pruning" << std::endl;

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

	prob->evaluate(box);
}

// TODO Replace this mock implementation with Acterberg's heuristic
void algorithm::lp_pruning() {

	// FIXME Check for convergence before tightening a variable
	for (int i=0; i<n; ++i) {

		box[i].tighten_bounds();
	}
}

void algorithm::delete_box() {

	delete[] box_orig;

	box_orig = 0;
}

void algorithm::rollback() {

	init_variables(box, box_orig, n);
}

void algorithm::split() {

	assert(false);
}

// TODO Replace mock implementation
bool algorithm::sufficient_progress() {

	copy_bounds(box, box_orig, n);

	return true;
}

}

