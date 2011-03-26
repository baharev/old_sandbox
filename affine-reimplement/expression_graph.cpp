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
#include <iostream> // FIXME Remove when ready
#include "expression_graph.hpp"
#include "box_generator.hpp"
#include "delete_struct.hpp"
#include "demangle.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"
#include "gap_info.hpp"
#include "interval.hpp"
#include "problem_data.hpp"
#include "sol_tracker.hpp"
#include "vector_dump.hpp"

using namespace std;

namespace asol {

template <typename T>
extern const std::vector<primitive<T>*> convert(const std::vector<primitive<builder>*>& v);

template <typename T>
expression_graph<T>::expression_graph(const problem_data* problem, const DoubleArray2D& solutions) :

v          (problem->peek_index()),
n_vars     (problem->number_of_variables()),
primitives (convert<T>(problem->get_primitives())),
constants  (problem->get_numeric_constants().begin(), problem->get_numeric_constants().end()),
initial_box(problem->get_initial_box()),
index_sets (problem->get_index_sets()),
constraints(problem->get_constraints()),
apriori_sols(solutions),
tracker    (0),
orig       (v.size()),
hull       (v.size())

{
	set_variables();

	primitive<T>::set_vector(&v);
}

template <typename T>
expression_graph<T>::~expression_graph() {

	for_each(primitives.begin(), primitives.end(), Delete());

	delete tracker;
}

template <typename T>
void expression_graph<T>::set_variables() {

	for (int i=0; i<n_vars; ++i) { // Required by affine ctor

		const Bounds& bound = initial_box.at(i);

		v.at(i) = T(bound.first, bound.second);
	}

	set_non_variables();

	set_numeric_consts();
}

template <typename T>
void expression_graph<T>::set_non_variables() {

	const double DMAX(1.0e+150);

	fill(v.begin()+n_vars, v.end(), T(-DMAX, DMAX));
}

template <typename T>
void expression_graph<T>::set_numeric_consts() {

	for_each(constants.begin(), constants.end(), bind1st(mem_fun(&expression_graph<T>::set_constant), this));
}

template <typename T>
void expression_graph<T>::set_constant(const Pair p) {

	ASSERT2(p.first>=n_vars, "index, n_vars: "<<p.first<<", "<<n_vars);

	v.at(p.first) = T(p.second);
}

template <typename T>
void expression_graph<T>::set_box(const T* box, const int length) {

	ASSERT(length == n_vars);

	copy(box, box+n_vars, v.begin()); // v.assign cannot be used here

	set_non_variables();

	// TODO Could save this if intersect would not create empty an interval
	set_numeric_consts();
}

template <typename T>
const T* expression_graph<T>::get_box() const {

	return &(v.at(0));
}

template <typename T>
void expression_graph<T>::evaluate_all() {

	for_each(primitives.begin(), primitives.end(), mem_fun(&primitive<T>::evaluate));
}

template <typename T>
void expression_graph<T>::revise_all() {

	evaluate_all();

	for_each(primitives.rbegin(), primitives.rend(), mem_fun(&primitive<T>::revise));
}

template <typename T>
int expression_graph<T>::constraint_begin(int i) const {

	return (i!=0) ? constraints.at(i-1)+1 : 0;
}

template <typename T>
int expression_graph<T>::constraint_end(int i) const {

	return constraints.at(i);
}

template <typename T>
void expression_graph<T>::evaluate_constraint(int k) {

	const int begin = constraint_begin(k);

	const int end   = constraint_end(k);

	for (int i=begin; i<=end; ++i) {

		primitives.at(i)->evaluate();

		//check_transitions_since_last_call();
	}
}

template <typename T>
void expression_graph<T>::revise_constraint(int k) {

	const int end   = constraint_end(k);

	const int begin = constraint_begin(k);

	for (int i=end; i>=begin; --i) {

		primitives.at(i)->revise();

		//check_transitions_since_last_call();
	}
}

template <typename T>
void expression_graph<T>::evaluate_up_to(const int k) {

	for (int i=0; i<=k; ++i) {

		evaluate_constraint(i);

		//check_transitions_since_last_call();
	}
}

template <typename T>
void expression_graph<T>::revise_up_to(const int k) {

	evaluate_up_to(k);

	for (int i=k; i>=0; --i) {

		revise_constraint(i);

		//check_transitions_since_last_call();
	}
}

template <typename T>
void expression_graph<T>::revise_all2() {

	const int last = constraints_size()-1;

	revise_up_to(last);
}

template <typename T>
void expression_graph<T>::iterative_revision() {

	const int end = constraints_size();

	for (int pos=0; pos<end; ++pos) {

		revise_up_to(pos);
	}
}

template <typename T>
struct raii {

	raii(vector<gap_info<T> >* gaps) {

		gaps->clear();

		primitive<T>::set_gap_container(gaps);
	}

	~raii() { primitive<T>::set_gap_container(0); }
};

template <typename T>
void expression_graph<T>::iterative_revision_save_gaps() {

	const int end = constraints_size();

	for (int pos=0; pos<end-1; ++pos) {

		revise_up_to(pos);
	}

	evaluate_up_to(end-1);

	raii<T> set_primitive_gap_container(&gaps);

	for_each(primitives.rbegin(), primitives.rend(), mem_fun(&primitive<T>::revise));

}

template <typename T>
void expression_graph<T>::probing() {

	revise_all(); // TODO Will it be always used like this?

	save_current_as_orig();

	const int m = constraints_size();

	for (int i=0; i<m; ++i) {

		hull.clear();

		set_orig_as_v();

		probe_in_constraint(i);
	}

	// the result of probing is in orig, result is retrieved as v
	set_orig_as_v();
}

template <typename T>
void expression_graph<T>::probe_in_constraint(const int k) {

	box_generator generator(v, index_sets.at(k), 3);

	if (generator.empty()) {
		// No progress, orig was set in the previous iteration
		return;
	}

	while (generator.get_next()) {

		set_orig_as_v();

		generator.set_box();

		revise_up_to_with_hull_saved(k);
	}

	write_hull_to_orig();
}

template <typename T>
void expression_graph<T>::save_current_as_orig() {

	orig.assign(v.begin(), v.end());
}

template <typename T>
void expression_graph<T>::set_orig_as_v() {

	v.assign(orig.begin(), orig.end());
}

template <typename T>
void expression_graph<T>::revise_up_to_with_hull_saved(const int k) {

	try {

		//save_containment_info();

		//show_variables(cout);

		revise_up_to(k);

		//show_variables(cout);

		//check_transitions_since_last_call();
	}
	catch (infeasible_problem& ) {
		// Print something?
		return;
	}
	catch (numerical_problems& ) {
		ASSERT2(false, "implementation not updated properly");
	}

	save_hull();
}

template <typename T>
void expression_graph<T>::save_hull() {

	if (hull.empty()) {

		hull.assign(v.begin(), v.end());
	}
	else {

		transform(hull.begin(), hull.end(), v.begin(), hull.begin(), hull_of);
	}
}

template <typename T>
void expression_graph<T>::write_hull_to_orig() {

	if (hull.empty()) {

		throw infeasible_problem();
	}

	bool changed = intersect_hull_and_orig();

	if (changed) {

		set_orig_as_v();

		iterative_revision_save_gaps();

		save_current_as_orig();
	}
}

template <typename T>
bool expression_graph<T>::intersect_hull_and_orig() {

	bool changed = false;

	try {

		changed = compute_intersection();
	}
	catch (infeasible_problem& ) {

		ASSERT(false);
	}

	return changed;
}

template <typename T>
bool expression_graph<T>::compute_intersection() {

	typename vector<T>::const_iterator end = hull.end();

	typename vector<T>::iterator result = orig.begin();

	bool changed = false;

	for (typename vector<T>::const_iterator i=hull.begin(); i!=end; ++i, ++result) {

		if (result->intersect(*i)) {

			changed = true;
		}
	}

	return changed;
}

template <typename T>
void expression_graph<T>::save_containment_info() {

	delete tracker;

	tracker = new sol_tracker(apriori_sols);

	tracker->save_containment_info(&v);
}

template <typename T>
void expression_graph<T>::print_containment_statistics() const {

	tracker->print_containment_statistics();
}

template <typename T>
bool expression_graph<T>::contains_solution() const {

	return tracker->contains_solution();
}

template <typename T>
void expression_graph<T>::increment_found_solution_counters() {

	tracker->increment_found_solution_counters();
}

template <typename T>
void expression_graph<T>::print_found_solutions() const {

	tracker->print_found_solutions();
}

template <typename T>
void expression_graph<T>::check_transitions_since_last_call() {

	tracker->check_transitions_since_last_call(&v);
}

template <typename T>
void expression_graph<T>::dump(const char* filename) const {

	asol::dump(v, filename);
}

template <typename T>
void expression_graph<T>::dump_trackers_previous() const {

	tracker->dump_previous_v();
}

template <typename T>
void expression_graph<T>::load_from_previous_dump() {

	const size_t size = v.size();

	asol::load(v);

	ASSERT2(size==v.size(),"size mismatch after loading, size: "<<size<<", "<<v.size());
}

template <typename T>
void expression_graph<T>::show_primitives_and_constraints(std::ostream& out) const {

	out << "Dumping primitives" << endl;

	const int n = static_cast<int> (primitives.size());

	for (int i=0; i<n; ++i) {
		out << i << ":\t" << name(typeid(*primitives.at(i))) << endl;
	}
	out << endl;

	out << "Dumping constraint indices" << endl;

	const int n_con = static_cast<int> (constraints.size());

	for (int i=0; i<n_con; ++i) {
		out << i << ":\t" << constraints.at(i) << endl;
	}
}

template <typename T>
void expression_graph<T>::show_variables(ostream& out) const {

	for (int i=0; i<n_vars; ++i) {

		out << i << ": " << v.at(i) << endl;
	}
}

template <typename T>
const T& expression_graph<T>::last_value() const {

	return v.at(v.size()-1);
}

template class expression_graph<interval>;

}
