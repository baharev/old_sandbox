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

#include <vector>
#include "operations.hpp"
#include "primitives.hpp"

namespace asol {

template <typename T>
class expression_graph : private operations {

public:

	typedef std::vector<primitive*> Vector;

	expression_graph(int number_of_variables, const Vector& p);

	void set_variables(const T box[], const int length);

	void evaluate_all();

	const T& last_value() const;

	void evaluate_primitive(int i);

	~expression_graph();

private:

	typedef Vector::iterator itr;

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	void clear_all();

	virtual void addition(int z, int x, int y)      { v.at(z)=v.at(x)+v.at(y); }
	virtual void substraction(int z, int x, int y)  { v.at(z)=v.at(x)-v.at(y); }
	virtual void multiplication(int z, int x, int y){ v.at(z)=v.at(x)*v.at(y); }
	virtual void division(int z, int x, int y)      { v.at(z)=v.at(x)/v.at(y); }

	std::vector<T> v;
	Vector primitives;
};

template <typename T>
expression_graph<T>::expression_graph(int number_of_variables, const Vector& p){

	v.resize(number_of_variables);

	primitives = Vector(p);
}

template <typename T>
expression_graph<T>::~expression_graph() {

	clear_all();
}

template <typename T>
void expression_graph<T>::clear_all() {

	v.clear();

	for (itr i = primitives.begin(); i!=primitives.end(); ++i) {

		delete *i;
	}

	primitives.clear();
}

template <typename T>
void expression_graph<T>::set_variables(const T box[], const int length) {

	for (int i=0; i<length; ++i) {

		v.at(i) = box[i];
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

};

#endif // EXPRESSION_GRAPH_HPP_
