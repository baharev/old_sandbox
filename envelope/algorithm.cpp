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

#include <assert.h>
#include "algorithm.hpp"
#include "envelope.hpp"
#include "exceptions.hpp"
#include "problem.hpp"

namespace asol {

algorithm::algorithm(const problem* const p)
	: n(p->size()), prob(p), box_orig(new interval[n]), box(new var[n])
{
	depth = boxes_processed = 0;

	pending.push_back(p->initial_box());
}

algorithm::~algorithm() {

	delete[] box_orig;

	delete[] box;

	var::release_all();
}

void algorithm::run() {

	do {

		iteration_step();

	} while (!pending.empty());
}

void algorithm::iteration_step() {

	try {

		contracting_step();
	}
	catch (infeasible_problem& ) {

		return;
	}
	catch (numerical_problems& ) {

		rollback();
	}

	// TODO Check for convergence

	if (sufficient_progress()) {

		push_front();

		return;
	}

	split();
}

void algorithm::contracting_step() {

	increment_counters();

	get_topmost_box();

	init_vars();

	evaluate();

	lp_pruning();
}

void algorithm::increment_counters() {

	// TODO Track and log depth info too
	++boxes_processed;
}

void algorithm::get_topmost_box() {

	interval* const box_current = pending.front();

	pending.pop_front();

	for (int i=0; i<n; ++i) {
		box_orig[i] = box_current[i];
	}

	delete[] box_current;
}

void algorithm::init_vars() {

	var::reset();
	// TODO Make convenience function
	for (int i=0; i<n; ++i) {
		box[i] = var(box_orig[i].inf(), box_orig[i].sup());
	}
}

void algorithm::evaluate() {

	prob->evaluate(box);
}

// TODO Replace this mock implementation with Acterberg's heuristic
void algorithm::lp_pruning() {

	// TODO It cannot become infeasible, if does then throw numerical error instead
	for (int i=0; i<n; ++i) {

		box[i].tighten_bounds();
	}
}

void algorithm::rollback() {

	assert(false);
}

void algorithm::push_front() {

	interval* const box_contracted = new interval[n];

	for (int i=0; i<n; ++i) {

		box_contracted[i] = box_orig[i];
	}

	pending.push_front(box_contracted);
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

