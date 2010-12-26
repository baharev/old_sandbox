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

#include <limits>
#include <vector>
#include <assert.h>
#include "operations.hpp"
#include "primitives.hpp"
#include "typedefs.hpp"

namespace asol {

template <typename T>
class expression_graph : private operations {

public:

	expression_graph(int number_of_variables,
			const PrimVector& primitives,
			const PairVector& numeric_constants);

	void set_variables(const T box[], const int length);

	void evaluate_all();

	const T& last_value() const;

	void evaluate_primitive(int i);

	~expression_graph();

private:

	typedef PrimVector::iterator itr;

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	int v_size() const { return static_cast<int> (v.size()); }
	int constants_size() const { return static_cast<int> (constants.size()); }
	void set_non_variables(const int length);
	void set_numeric_consts(const int length);

	virtual void addition(int z, int x, int y);
	virtual void substraction(int z, int x, int y);
	virtual void multiplication(int z, int x, int y);
	virtual void division(int z, int x, int y);
	virtual void square(int z, int x);

	std::vector<T> v;
	PrimVector primitives;
	const PairVector constants;
};

template <typename T>
expression_graph<T>::expression_graph(int number_of_variables,
									  const PrimVector& p,
									  const PairVector& numeric_const)
: v(number_of_variables), primitives(p), constants(numeric_const)
{

}

template <typename T>
expression_graph<T>::~expression_graph() {

	for (itr i = primitives.begin(); i!=primitives.end(); ++i) {

		delete *i;
	}
}

template <typename T>
void expression_graph<T>::set_variables(const T box[], const int length) {

	for (int i=0; i<length; ++i) {

		v.at(i) = box[i];
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

		assert(index>=length);

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

};

#endif // EXPRESSION_GRAPH_HPP_
