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
#include "affine.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"
#include "expression_graph.hpp"
#include "index_recorder.hpp"
#include "splitting_strategy.hpp"
#include "interval.hpp"
#include "lp_solver.hpp"
#include "problem.hpp"
#include "problem_data.hpp"
#include "vector_dump.hpp"

using std::cout;
using std::endl;
using std::fabs;

namespace asol {

search_procedure::search_procedure(const problem<builder>* p)
: prob(p),
  //split_strategy(new max_diam_selector(prob->number_of_variables())),
  split_strategy(new Jacobsen_x1_D(prob->number_of_variables())),
  //split_strategy(new eco9_sparsity(prob->number_of_variables())),
  n_vars(prob->number_of_variables()),
  representation(0),
  lp(new lp_solver),
  box_orig(0)
{
	build_problem_representation();

	init_dags();

	init_lp_solver();

	solutions_found = splits = boxes_processed = 0;

	representation = 0;

	delete prob;

	prob = 0;

	builder::reset();
}

search_procedure::~search_procedure() {

	delete ia_dag;

	delete aa_dag;

	delete lp;

	delete split_strategy;

	lp_solver::free_environment();
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

void search_procedure::init_dags() {

	//const IntArray2D index_set = index_sets();

	ia_dag = new expression_graph<interval>(representation, prob->solutions());

	push_initial_box_to_deque();
	//dbg_initial_box_from_dump();

	affine::set_vector(ia_dag->get_v());

	aa_dag = new expression_graph<affine>(representation);
}

void search_procedure::init_lp_solver() {

	index_recorder rec(representation);

	lp->set_pruning_indices(rec.lp_pruning_index_sets());

	lp->set_number_of_vars(n_vars);

	lp->set_affine_vars(aa_dag->get_v());

	affine::set_lp_solver(lp);
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

std::vector<std::vector<int> > search_procedure::index_sets() const {

	//index_recorder rec(representation);

	//rec.dump();

	//return rec.constraint_index_sets();

	// FIXME Just a hack for probing2
	IntArray2D indices;

	IntVector vec;
	vec.push_back(0);
	vec.push_back(15);

	indices.push_back(vec);

	return indices;
}

void search_procedure::dbg_initial_box_from_dump() {

	ASSERT(pending_boxes.empty());

	std::vector<interval> v;

	load(v);

	interval* x = new interval[n_vars];

	std::copy(&v.at(0), &v.at(n_vars), x);

	pending_boxes.push_back(x);
}

void search_procedure::run() {

	while (has_more_boxes()) {

		get_next_box();

		process_box();

		ia_dag->show_variables(cout);

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

	lp->show_iteration_count();
	ia_dag->print_found_solutions();
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

		dbg_check_infeasibilty();
		delete_box();
	}
	catch (numerical_problems& ) {

		roll_back();
	}
	catch (convergence_reached& ) {

		print_box();
		dbg_solution_count();
		delete_box();
	}
}

void search_procedure::dbg_check_infeasibilty() const {

	if(ia_dag->contains_solution()) {

		ia_dag->dump_trackers_previous();

		ia_dag->dump("v_current_dump.txt");

		ASSERT(false);
	}
}

void search_procedure::dbg_solution_count() {

	ia_dag->increment_found_solution_counters();
}

void search_procedure::roll_back() {

	cout << "Warning: numerical problems, rolling back!" << endl;

	ia_dag->set_box(box_orig, n_vars);
}

void search_procedure::contracting_step() {

	ia_dag->set_box(box_orig, n_vars);

	ia_dag->save_containment_info();
	// TODO Check index sets!
	//ia_dag->probing2();

	ia_dag->iterative_revision();

	ia_dag->check_transitions_since_last_call();

	check_convergence();

	lp->reset();

	aa_dag->reset_vars();

	aa_dag->evaluate_all();

	lp->check_feasibility(); // FIXME Once found feasible, cannot become infeas!!!

	lp->prune(std::vector<int>()); // Cannot throw infeasible problem

	check_convergence();

	ia_dag->check_transitions_since_last_call();

	ia_dag->iterative_revision();

	ia_dag->check_transitions_since_last_call();

	check_convergence();
}

const double CONVERGENCE_TOL = 0.05; // FIXME Just for testing

struct wide {

	bool operator()(const interval& x) const { return !x.is_narrow(CONVERGENCE_TOL); }
};

void search_procedure::check_convergence() {

	// TODO Move convergence check to expression_graph?
	const interval* const box = ia_dag->get_box();

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

	ia_dag->show_variables(cout);
}

bool search_procedure::sufficient(const double max_progress) const {

	return max_progress < 0.95; // TODO Introduce options class
	                            // FIXME See also affine::excellent_progress_made
}

bool search_procedure::sufficient_progress() {

	cout << "Computing progress" << endl;

	const double best_reduction = compute_max_progress();

	if (sufficient(best_reduction)) {

		cout << "Sufficient progress made" << endl;
		cout << "-----------------------------------------------------" << endl;
	}

	const interval* const box = ia_dag->get_box();

	std::copy(box, box+n_vars, box_orig); // FIXME Sort of hidden side-effect

	return sufficient(best_reduction);
}

struct diam_reduction {

	double operator()(const interval& x, const interval& y) {

		ASSERT(x.subset_of(y));

		double result = 10; // TODO Magic number

		if (!x.is_narrow(CONVERGENCE_TOL) && y.diameter()!=0) {

			result = x.diameter() / y.diameter();
		}

		return result;
	}
};

double search_procedure::compute_max_progress() const {

	const interval* const box_contracted = ia_dag->get_box();

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

	const int index = split_strategy->index_to_split(box_orig);

	ASSERT2(0<=index && index < n_vars, "index: "<<index);

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

}
