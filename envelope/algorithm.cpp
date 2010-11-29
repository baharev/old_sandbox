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

#include "algorithm.hpp"
#include "envelope.hpp"
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

	get_topmost_box();

	do {

		init_vars();

		build_lp();

		// TODO Achterberg's contraction

	} while (sufficient_progress());
}

void algorithm::get_topmost_box() {

	interval* const box_current = pending.front();

	pending.pop_front();

	for (int i=0; i<n; ++i) {
		box_orig[i] = box_current[i];
	}
}

void algorithm::init_vars() {

	var::reset();
	// TODO Make convenience function
	for (int i=0; i<n; ++i) {
		box[i] = var(box_orig[i].inf(), box_orig[i].sup());
	}
}

void algorithm::build_lp() {

	prob->evaluate(box);
}

// TODO Replace this mock implementation
bool algorithm::sufficient_progress() {

	// TODO Implement convenience function to copy new bounds

	for (int i=0; i<n; ++i) {

		box_orig[i] = box[i].compute_bounds();
	}

	return true;
}

}

