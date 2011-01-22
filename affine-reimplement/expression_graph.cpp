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
#include <limits>
#include "expression_graph.hpp"
#include "delete_struct.hpp"
#include "diagnostics.hpp"
#include "interval.hpp"
#include "problem_data.hpp"

using namespace std;

namespace asol {

template <typename T>
extern const std::vector<primitive<T>*> convert(const std::vector<primitive<builder>*>& v);

template <typename T>
expression_graph<T>::expression_graph(const problem_data* problem) :

v          (problem->peek_index()),
primitives (convert<T>(problem->get_primitives())),
constants  (problem->get_numeric_constants()),
initial_box(problem->get_initial_box()),
index_sets (problem->get_index_sets()),
constraints(problem->get_constraints())

{
	set_variables();
	primitive<T>::set_vector(&v);
}

template <typename T>
expression_graph<T>::~expression_graph() {

	for_each(primitives.begin(), primitives.end(), Delete());
}

template <typename T>
void expression_graph<T>::set_variables() {

	const int length = static_cast<int> (initial_box.size());

	for (int i=0; i<length; ++i) {

		const Bounds& bound = initial_box.at(i);

		v.at(i) = T(bound.first, bound.second);
	}

	set_non_variables(length);

	set_numeric_consts(length);
}

template <typename T>
void expression_graph<T>::set_non_variables(const int length) {

	const double DMAX =  numeric_limits<double>::max();
	const double DMIN = -DMAX;

	for (int i=length; i<v_size(); ++i) {

		v.at(i) = T(DMIN, DMAX);
	}
}

template <typename T>
void expression_graph<T>::set_numeric_consts(const int length) {

	Map::const_iterator i = constants.begin();

	while (i!=constants.end()) {

		const int    index = i->first;
		const double value = i->second;

		ASSERT2(index>=length, "index, length: "<<index<<", "<<length);

		v.at(index) = T(value);

		++i;
	}
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

	for (int i=k; i>=0; --i) {

		revise_constraint(i);
	}
}

template <typename T>
void expression_graph<T>::revise_all2() {

	const int last = constraints_size()-1;

	evaluate_up_to(last);

	revise_up_to(last);
}

template <typename T>
void expression_graph<T>::probing() {

	const int end = constraints_size();

	for (int pos=0; pos<end; ++pos) {

		evaluate_up_to(pos);

		revise_up_to(pos);
	}
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

template <typename T>
void expression_graph<T>::evaluate_primitive(int i) {

	primitives.at(i)->evaluate();
}

template class expression_graph<interval>;

}
