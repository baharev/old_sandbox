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
#include "expression_graph.hpp"
#include "box_generator.hpp"
#include "delete_struct.hpp"
#include "diagnostics.hpp"
#include "exceptions.hpp"
#include "interval.hpp"
#include "problem_data.hpp"

using namespace std;

namespace asol {

template <typename T>
extern const std::vector<primitive<T>*> convert(const std::vector<primitive<builder>*>& v);

template <typename T>
expression_graph<T>::expression_graph(const problem_data* problem) :

v          (problem->peek_index()),
n_vars     (problem->number_of_variables()),
primitives (convert<T>(problem->get_primitives())),
constants  (problem->get_numeric_constants().begin(), problem->get_numeric_constants().end()), // TODO Do it in problem_data
initial_box(problem->get_initial_box()),
index_sets (problem->get_index_sets()),
constraints(problem->get_constraints())

{
	set_variables();
	primitive<T>::set_vector(&v);

	orig.resize(v.size());
	hull.resize(v.size());
}

template <typename T>
expression_graph<T>::~expression_graph() {

	for_each(primitives.begin(), primitives.end(), Delete());
}

template <typename T>
void expression_graph<T>::set_variables() {

	for (int i=0; i<n_vars; ++i) {

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
	}
}

template <typename T>
void expression_graph<T>::revise_constraint(int k) {

	const int end   = constraint_end(k);

	const int begin = constraint_begin(k);

	for (int i=end; i>=begin; --i) {

		primitives.at(i)->revise();
	}
}

template <typename T>
void expression_graph<T>::evaluate_up_to(const int k) {

	for (int i=0; i<=k; ++i) {

		evaluate_constraint(i);
	}
}

template <typename T>
void expression_graph<T>::revise_up_to(const int k) {

	evaluate_up_to(k);

	for (int i=k; i>=0; --i) {

		revise_constraint(i);
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
void expression_graph<T>::probing() {

	save_current_as_orig();

	const int m = constraints_size();

	for (int i=0; i<m; ++i) {

		probing(i);
	}
}

template <typename T>
void expression_graph<T>::probing(const int k) {

	hull.clear();

	set_orig_as_v();

	box_generator generator(v, index_sets.at(k), 3);

	if (generator.empty()) {

		return;
	}

	while (generator.get_next()) {

		set_orig_as_v();

		generator.set_box();

		probe_one(k);
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
void expression_graph<T>::probe_one(const int k) {

	try {

		revise_up_to(k);
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

		revise_all();

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
void expression_graph<T>::show_variables(ostream& out) const {

	const int length = static_cast<int> (initial_box.size());

	for (int i=0; i<length; ++i) {

		out << i << ": " << v.at(i) << endl;
	}
}

template <typename T>
const T& expression_graph<T>::last_value() const {

	return v.at(v.size()-1);
}

template class expression_graph<interval>;

}
