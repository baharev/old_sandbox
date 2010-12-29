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

#include <iosfwd>
#include <limits>
#include <vector>
#include "operations.hpp"
#include "diagnostics.hpp"
#include "primitives.hpp"
#include "typedefs.hpp"

namespace asol {

template <typename T>
class expression_graph : private operations {

public:

	expression_graph(int number_of_arguments,
			const PrimVector& primitives,
			const PairVector& numeric_constants,
			const PairVector& constraint_rhs,
			const BoundVector& initialbox);

	void evaluate_all();

	void revise_all();

	void show_variables(std::ostream& out) const;

	const T& last_value() const;

	void evaluate_primitive(int i);

	~expression_graph();

private:

	typedef PrimVector::iterator itr;
	typedef PrimVector::reverse_iterator ritr;

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	int v_size() const { return static_cast<int> (v.size()); }
	int constants_size() const { return static_cast<int> (constants.size()); }
	void set_variables();
	void set_non_variables(const int length);
	void set_numeric_consts(const int length);

	virtual void addition(int z, int x, int y);
	virtual void substraction(int z, int x, int y);
	virtual void multiplication(int z, int x, int y);
	virtual void division(int z, int x, int y);
	virtual void square(int z, int x);
	virtual void exponential(int z, int x);
	virtual void equality_constraint(int body, int rhs);

	virtual void addition_revise(int z, int x, int y);
	virtual void substraction_revise(int z, int x, int y);
	virtual void multiplication_revise(int z, int x, int y);
	virtual void division_revise(int z, int x, int y);
	virtual void square_revise(int z, int x);
	virtual void exponential_revise(int z, int x);
	virtual void equality_constraint_revise(int body, int rhs);

	std::vector<T> v;
	PrimVector primitives;
	const PairVector constants;
	const PairVector rhs_constraints;
	const BoundVector initial_box;
};

template <typename T>
expression_graph<T>::expression_graph(int number_of_arguments,
									  const PrimVector& p,
									  const PairVector& numeric_const,
									  const PairVector& constraint_rhs,
									  const BoundVector& initialbox)
: v(number_of_arguments),
  primitives(p),
  constants(numeric_const),
  rhs_constraints(constraint_rhs),
  initial_box(initialbox)
{
	set_variables();
}

template <typename T>
expression_graph<T>::~expression_graph() {

	for (itr i = primitives.begin(); i!=primitives.end(); ++i) {

		delete *i;
	}
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

	for (int i=0; i<constants_size(); ++i) {

		const int    index = constants.at(i).first;
		const double value = constants.at(i).second;

		ASSERT2(index>=length, "index, length: "<<index<<", "<<length)

		v.at(index) = T(value);
	}
}

template <typename T>
void expression_graph<T>::evaluate_all() {

	for (itr i=primitives.begin(); i!=primitives.end(); ++i) {

		(*i)->evaluate(this);
	}
}

template <typename T>
void expression_graph<T>::revise_all() {

	evaluate_all();

	for (ritr i=primitives.rbegin(); i!=primitives.rend(); ++i) {

		(*i)->revise(this);
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

	primitives.at(i)->evaluate(this);
}

template <typename T>
void expression_graph<T>::addition(int z, int x, int y) {

	v.at(z).assign( v.at(x)+v.at(y) );
}

template <typename T>
void expression_graph<T>::substraction(int z, int x, int y) {

	v.at(z).assign( v.at(x)-v.at(y) );
}

template <typename T>
void expression_graph<T>::multiplication(int z, int x, int y) {

	v.at(z).assign( v.at(x)*v.at(y) );
}

template <typename T>
void expression_graph<T>::division(int z, int x, int y) {

	v.at(z).assign( v.at(x)/v.at(y) );
}

template <typename T>
void expression_graph<T>::square(int z, int x) {

	v.at(z).assign( sqr(v.at(x)) );
}

template <typename T>
void expression_graph<T>::exponential(int z, int x) {

	v.at(z).assign( exp(v.at(x)) );
}

template <typename T>
void expression_graph<T>::equality_constraint(int body, int rhs) {

	const double rhs_value = rhs_constraints.at(rhs).second;

	v.at(body).equals(rhs_value);
}

template <typename T>
void expression_graph<T>::addition_revise(int z, int x, int y) {

	addition_inverse(v.at(z), v.at(x), v.at(y));
}

template <typename T>
void expression_graph<T>::substraction_revise(int z, int x, int y) {

	substraction_inverse(v.at(z), v.at(x), v.at(y));
}

template <typename T>
void expression_graph<T>::multiplication_revise(int z, int x, int y) {

	multiplication_inverse(v.at(z), v.at(x), v.at(y));
}

template <typename T>
void expression_graph<T>::division_revise(int z, int x, int y) {

	division_inverse(v.at(z), v.at(x), v.at(y));
}

template <typename T>
void expression_graph<T>::square_revise(int z, int x) {

	sqr_inverse(v.at(z), v.at(x));
}

template <typename T>
void expression_graph<T>::exponential_revise(int z, int x) {

	exp_inverse(v.at(z), v.at(x));
}

// TODO Not yet clear how to invert a constraint
template <typename T>
void expression_graph<T>::equality_constraint_revise(int body, int rhs) {

	const double rhs_value = rhs_constraints.at(rhs).second;

	T rhs_val(rhs_value);

	equality_constraint_inverse(v.at(body), rhs_val);
}

};

#endif // EXPRESSION_GRAPH_HPP_
