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
//#include <iomanip>
#include <cmath> // FIXME
#include <assert.h>
#include "algorithm.hpp"
#include "constants.hpp"
#include "envelope.hpp"
#include "exceptions.hpp"
#include "problem.hpp"

using std::cout;
using std::endl;

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

		set_topmost_box();

		iteration_step();

	} while (has_more_boxes());

	print_statistics();
}

bool algorithm::has_more_boxes() const {

	return !pending.empty();
}

void algorithm::add_initial_box() {

	pending.push_back(prob->initial_box());
}

void algorithm::set_topmost_box() {

	cout << "=========================================================" << endl;

	assert(box_orig==0);

	box_orig = pending.front();

	pending.pop_front();
}

void algorithm::iteration_step() {

	bool deleted = false;

	do {

		deleted = one_pass();

	} while ( !deleted && sufficient_progress() );

	if(!deleted) {

		split();
	}
}

bool algorithm::one_pass() {

	bool deleted = false;

	try {
		contracting_step();
	}
	catch (infeasible_problem& ) {
		delete_box();
		deleted = true;
	}
	catch (numerical_problems& ) {
		rollback();
	}
	catch (convergence_reached& ) {
		print_box();
		delete_box();
		deleted = true;
	}

	return deleted;
}

void algorithm::contracting_step() {

	increment_counters();

	initialize_variables();

	evaluate();

	try {

		lp_pruning();
	}
	catch (infeasible_problem& ) {

		cout << "Warning: numerical problems in LP pruning" << endl;

		throw numerical_problems();
	}
}

void algorithm::increment_counters() {

	// TODO Track and log depth info too
	++boxes_processed;
}

void algorithm::initialize_variables() {

	init_variables(box, box_orig, n);
}

void algorithm::evaluate() {

	cout << "Evaluation" << endl;

	prob->evaluate(box);

	check_convergence();
}

void algorithm::lp_pruning() {

	cout << "Running LP pruning" << endl;

	prune_all(box, n);

	check_convergence();
}

void algorithm::delete_box() {

	cout << "Box discarded" << endl;

	delete[] box_orig;

	box_orig = 0;
}

void algorithm::rollback() {

	cout << "Warning: numerical problems detected" << endl;

	init_variables(box, box_orig, n);
}

void algorithm::print_box() const {

	for (int i=0; i<n; ++i) {

		cout << box[i] << endl;
	}
}

void algorithm::check_convergence() {

	int index = find_max_width(box, n);

	if (box[index].width() <= TOL_SOLVED) {

		cout << "Found a solution!" << endl;

		++solutions_found;

		throw convergence_reached();
	}
}
/*
void algorithm::dbg_box_width() const {

	using namespace std;

	bool to_dump = false;

	for (int i=0; i<n; ++i) {
		double width1 = box[i].width();
		double width2 = box_orig[i].diameter();
		if (width1!=width2) {
			width1 = width2 = 0;
			to_dump = true;
			break;
		}
	}

	if (to_dump) {
		cout << "###" << endl;
		for (int i=0; i<n; ++i) {
			cout << setprecision(16) << scientific << box[i] << endl;
			cout << setprecision(16) << scientific << box_orig[i] << endl;
			cout << i << ": " << (box[i].width()==box_orig[i].diameter()) << endl;
			bool lb_equals = box[i].range.lb == box_orig[i].lb;
			bool ub_equals = box[i].range.ub == box_orig[i].ub;
			cout << lb_equals << "  " << ub_equals << endl;
		}
	}
}
*/
int algorithm::select_index_to_split() const {

	//double x1 = box[0].width();
	//double D  = box[15].width();

	//int index = (x1 > D)? 0 : 15;

	int index = find_max_width(box, n);

	cout << "Splitting " << index << ", " << box_orig[index] << endl;
	cout << "var[index] = " << box[index] << endl;

	// FIXME They may not equal exactly due to inlining?
	assert ( std::fabs(box[index].width()- box_orig[index].diameter()) < 1.0e-6);
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

double algorithm::compute_max_progress(const interval box_contracted[]) const {

	double best_reduction = 10;

	for (int i=0; i<n; ++i) {

		if (box_orig[i].diameter() == 0) {
			cout << "Warning: index " << i << "has zero width" << endl;
			continue;
		}

		double reduction = box_contracted[i].diameter() / box_orig[i].diameter();

		if (reduction<best_reduction) {
			best_reduction = reduction;
		}
	}

	cout << "Best reduction: " << best_reduction << endl;

	return best_reduction;
}

bool algorithm::sufficient_progress() {

	cout << "Computing progress" << endl;

	interval box_contracted[n];

	copy_bounds(box, box_contracted, n);

	double best_reduction = compute_max_progress(box_contracted);

	bool sufficient = false;

	if (best_reduction < 0.95) {

		sufficient = true;

		cout << "Sufficient progress made" << endl;
		cout << "-----------------------------------------------------" << endl;
	}

	copy_bounds(box, box_orig, n);

	//dbg_box_width();

	return sufficient;
}

void algorithm::print_statistics() const {

	cout << endl;
	cout << "=========================================================" << endl;
	cout << "Number of splits: " << splits << ", solutions: ";
	cout << solutions_found << endl;
}

}

