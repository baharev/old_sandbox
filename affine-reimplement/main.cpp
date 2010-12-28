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
#include <exception>
#include <typeinfo>
#include "builder.hpp"
#include "diagnostics.hpp"
#include "expression_graph.hpp"
#include "interval.hpp"
#include "problem.hpp"

using namespace asol;

//==============================================================================

template <typename T>
class Hansen_example : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum { X, Y, SIZE };

};

template <typename T>
int Hansen_example<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Hansen_example<T>::initial_box() const {

	T* box = new T[SIZE];

	//box[X] = T(1.33073);
	//box[Y] = T(1);

	//box[X] = T(1);
	//box[Y] = T(10);

	box[X] = T(1, 10);
	box[Y] = T(1, 10);

	return box;
}

template <typename T>
void Hansen_example<T>::evaluate(const T v[]) const {

	const T& x = v[X];
	const T& y = v[Y];

	T xy = x*y;

	// FIXME Expression graph does not call common subexpression on interval
	xy.mark_as_common_subexpression();

	T z = (5*x-4*sqr(y)+14*xy)/(sqr(x)+y+xy);

	z.equals(10);
}

//==============================================================================

template <typename T>
class Example_challange : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum { W, X, Y, Z, B, C, A, SIZE };

};

template <typename T>
int Example_challange<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Example_challange<T>::initial_box() const {

	T* box = new T[SIZE];

	box[W] = T(-0.9, -0.6);
	box[X] = T(-0.1,  0.2);

	box[Y] = T( 0.3,  0.7);
	box[Z] = T(-0.2,  0.1);

	box[B] = T(-1, 1);
	box[C] = T(-1, 1);

	box[A] = T(7.0, 9.0);

	return box;
}

template <typename T>
void Example_challange<T>::evaluate(const T var[]) const {

	const T& w = var[W];
	const T& x = var[X];

	const T& y = var[Y];
	const T& z = var[Z];

	const T& b = var[B];
	const T& c = var[C];

	const T& a = var[A];

	const T u = sqr(w) + sqr(x);

	u.mark_as_common_subexpression();

	const T v = sqr(y) + sqr(z);

	v.mark_as_common_subexpression();

	const T d = b*(x*y - w*z) + c*(x*z + w*y);

	const T f = (a*(u-v)+2*d)/(u+v);

}

//==============================================================================

template <typename T>
class Example_1 : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum { X, SIZE };

};

template <typename T>
int Example_1<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Example_1<T>::initial_box() const {

	T* box = new T[SIZE];

	box[X] = T(2.0, 4.0);

	return box;
}

template <typename T>
void Example_1<T>::evaluate(const T v[]) const {

	const T& x = v[X];

	const T  y = (x-1)/(sqr(x)+2);
}

//==============================================================================

template <typename T>
class Example_2 : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum { X, SIZE };

};

template <typename T>
int Example_2<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Example_2<T>::initial_box() const {

	T* box = new T[SIZE];

	box[X] = T(1.0, 2.0);

	return box;
}

template <typename T>
void Example_2<T>::evaluate(const T v[]) const {

	const T& x = v[X];

	const T  y = (sqr(x)+x)/(16*x-9);
}

//==============================================================================

void build(const problem<builder>* prob) {

	builder::reset();

	builder* box = prob->initial_box();

	prob->evaluate(box);

	delete[] box;

	ASSERT(builder::number_of_variables() == prob->number_of_variables())

	delete prob;

	builder::record_occurence_info();
}

void dag_test(const problem<builder>* prob) {

	build(prob);

	expression_graph<interval> dag( builder::number_of_arguments(),
									builder::get_primitives(),
									builder::get_numeric_constants(),
									builder::get_rhs_of_constraints(),
									builder::get_initial_box());

	builder::dbg_dump_type_of_primitives();
	builder::dbg_show_info();
	builder::reset();

	dag.evaluate_all();

	std::cout << "Last value: " << dag.last_value() << std::endl;
}

void example_Hansen() {

	dag_test(new Hansen_example<builder> ());
}

void example_challenge() {

	dag_test(new Example_challange<builder> ());
}


void example_1() {

	dag_test(new Example_1<builder> ());
}


void example_2() {

	dag_test(new Example_2<builder> ());
}

void assert_tests() {

	double lb(0), ub(2), x(1);

	ASSERT2(lb<=x && x<=ub, "lb, x, ub: "<<lb<<", "<<x<<", "<<ub)

	ASSERT(lb<=x && x<=ub)

	interval a(1, 3), b(-2, 5);

	interval c;

	//c.diameter();

	a/b;
}

int main() {

	try {

		assert_tests();
	}
	catch (std::exception& e) {

		std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
	}

	example_Hansen();

	example_challenge();

	example_1();

	example_2();

	return 0;
}
