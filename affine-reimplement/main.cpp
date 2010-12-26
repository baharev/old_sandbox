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

#include <iostream>
#include "builder.hpp"
#include "expression_graph.hpp"
#include "interval.hpp"

using namespace asol;

enum { X, Y, SIZE };

int size() {

	return SIZE;
}

template <typename T>
T* initial_box() {

	T* box = new T[SIZE];

	box[X] = T(1.33073);
	box[Y] = T(1);

	return box;
}

template <typename T>
void f(const T v[]) {

	const T& x = v[X];
	const T& y = v[Y];

	T xy = x*y;

	T z = (5*x-4*sqr(y)+14*xy)/(sqr(x)+y+xy);
}

void build() {

	builder* box = initial_box<builder>();

	f(box);

	delete[] box;
}

template <typename T>
void set_box(expression_graph<T>& dag) {

	T* box = initial_box<T>();

	dag.set_variables(box, size());

	delete[] box;
}

void Main() {

	build();

	expression_graph<interval> dag(builder::number_of_variables(), builder::get_primitives());

	builder::reset();

	set_box(dag);

	dag.evaluate_all();

	std::cout << "Last value: " << dag.last_value() << std::endl;
}

int main() {

	Main();

	return 0;
}
