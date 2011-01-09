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

#ifndef EXPRESSION_GRAPH_HPP_
#define EXPRESSION_GRAPH_HPP_

#include <algorithm>
#include <iosfwd>
#include <limits>
#include <map>
#include <vector>
#include "diagnostics.hpp"
#include "primitives.hpp"
#include "typedefs.hpp"

namespace asol {

template <typename T>
class expression_graph {

public:

	typedef std::vector<primitive<T>*> PrimVector;

	expression_graph(int number_of_arguments,
			const PrimVector& primitives,
			const std::map<int,double>& numeric_constants,
			const PairVector& constraint_rhs,
			const BoundVector& initialbox);

	void evaluate_all();

	void revise_all();

	void show_variables(std::ostream& out) const;

	const T& last_value() const;

	void evaluate_primitive(int i);

	~expression_graph();

private:

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	int v_size() const { return static_cast<int> (v.size()); }
	void set_variables();
	void set_non_variables(const int length);
	void set_numeric_consts(const int length);

	std::vector<T> v;
	const PrimVector primitives;
	const std::map<int,double> constants;
	const PairVector rhs_constraints;
	const BoundVector initial_box;
};

template <typename T>
expression_graph<T>::expression_graph(int number_of_arguments,
									  const PrimVector& p,
									  const std::map<int,double>& numeric_const,
									  const PairVector& constraint_rhs,
									  const BoundVector& initialbox)
: v(number_of_arguments),
  primitives(p),
  constants(numeric_const),
  rhs_constraints(constraint_rhs),
  initial_box(initialbox)
{
	set_variables();
	primitive<T>::set_vector(&v);
}

template <typename T>
expression_graph<T>::~expression_graph() {

	std::for_each(primitives.begin(), primitives.end(), Delete());
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

	const double DMAX =  std::numeric_limits<double>::max();
	const double DMIN = -DMAX;

	for (int i=length; i<v_size(); ++i) {

		v.at(i) = T(DMIN, DMAX);
	}
}

template <typename T>
void expression_graph<T>::set_numeric_consts(const int length) {

	std::map<int,double>::const_iterator i = constants.begin();

	while (i!=constants.end()) {

		const int    index = i->first;
		const double value = i->second;

		ASSERT2(index>=length, "index, length: "<<index<<", "<<length)

		v.at(index) = T(value);

		++i;
	}
}

template <typename T>
void expression_graph<T>::evaluate_all() {
	// TODO Replace with for_each
	typename PrimVector::const_iterator i = primitives.begin();

	for (; i!=primitives.end(); ++i) {

		(*i)->evaluate();
	}
}

template <typename T>
void expression_graph<T>::revise_all() {

	evaluate_all();

	typename PrimVector::const_reverse_iterator i = primitives.rbegin();

	for (; i!=primitives.rend(); ++i) {

		(*i)->revise();
	}
}

template <typename T>
void expression_graph<T>::show_variables(std::ostream& out) const {

	const int length = static_cast<int> (initial_box.size());

	for (int i=0; i<length; ++i) {

		out << i << ": " << v.at(i) << std::endl;
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

};

#endif // EXPRESSION_GRAPH_HPP_
