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
#include "demangle.hpp"
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

	const T xy = x*y;

	// FIXME Expression graph does not call common subexpression on interval
	xy.mark_as_common_subexpression();

	const T z = (5*x-4*sqr(y)+14*xy)/(sqr(x)+y+xy);

	z.equals(14.5);
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

	f.equals(0);
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

template <typename T>
class Example_3 : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum { X, Y, SIZE };

};

template <typename T>
int Example_3<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Example_3<T>::initial_box() const {

	T* box = new T[SIZE];

	//box[X] = T( 0.0, 1.0);
	//box[Y] = T(-1.0, 1.0);

	box[X] = T(-1.0,-0.5);
	box[Y] = T( 0.5, 1.0);

	return box;
}

template <typename T>
void Example_3<T>::evaluate(const T v[]) const {

	const T& x = v[X];
	const T& y = v[Y];

	const T x2 = sqr(x);

	x2.mark_as_common_subexpression();

	const T eq1 = x2 + sqr(y);

	eq1.equals(1);

	const T eq2 = x2 - y;

	eq2.equals(0);
}

//==============================================================================

template <typename T>
class Jacobsen : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum {
		X1, X2, X3, X4, X5, X6, X7, X8,
		v1, v2, v3, v4, v5, v6, v7, C, SIZE
	};

};

template <typename T>
int Jacobsen<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Jacobsen<T>::initial_box() const {

	T* v = new T[SIZE];

	//      T(1.0e-4, 1.0);
	v[X1] = T(0.93, 0.94);

	for (int i=X2; i<=X8; ++i) {
		v[i] = T(1.0e-4, 1.0);
	}

	for (int i=v1; i<=v7; ++i) {
		v[i] = T(2.0, 4.0);
	}
	//     T(0.0, 1.12);
	v[C] = T(0.50, 0.51);

	return v;
}

// TODO Make an interval function for this: it is monotonous!
template <typename T>
const T H_Vap(const T& x) {

	return 0.1349*exp(-3.98*x) + 0.4397*exp(-0.088*x);
}

template <typename T>
const T H_Liq(const T& x) {

	return 0.1667*exp(-1.087*x);
}

template <typename T>
const T y_eq(const T& x) {

	return 3.55/( 1/x + 2.55 );
}

template <typename T>
void Jacobsen<T>::evaluate(const T v[]) const {

	const T& x1 = v[X1];
	const T& x2 = v[X2];
	const T& x3 = v[X3];
	const T& x4 = v[X4];
	const T& x5 = v[X5];
	const T& x6 = v[X6];
	const T& x7 = v[X7];
	const T& x8 = v[X8];
	const T& V1 = v[v1];
	const T& V2 = v[v2];
	const T& V3 = v[v3];
	const T& V4 = v[v4];
	const T& V5 = v[v5];
	const T& V6 = v[v6];
	const T& V7 = v[v7];
	const T& D = v[C];

	const T y1 = y_eq(x1);

	const T d = D*y1;

	d.mark_as_common_subexpression();

	const T Lw = (V1 - D)*(0.6010 - 0.2806*y1);

	Lw.equals(0.96);

	//--------------------------------------------------------------------------

	const T HV1 = H_Vap(x1);

	const T HL0 = H_Liq(y1);

	const T Q = V1*(HV1 - HL0) + D*HL0;

	Q.mark_as_common_subexpression();

	//==========================================================================

	const T M8 = (1.0-D)*x8 + d;

	M8.equals(0.5);

	//==========================================================================

	const T L7 = 4.0-D;

	const T y8 = y_eq(x8);

	const T M7 = 3.0*y8-L7*x7-d;

	M7.equals(-0.5);

	//--------------------------------------------------------------------------

	const T HV8 = H_Vap(x8);

	const T HL7 = H_Liq(x7);

	const T H7 = 3*HV8 -L7*HL7 - Q;

	H7.equals(-0.0968047);

	//==========================================================================

	const T y2 = y_eq(x2);

	const T M1 = V2* y2 - (    V2 - D)* x1 - d;

	M1.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV2 = H_Vap(x2);

	const T HL1 = H_Liq(x1);

	const T H1 = V2*HV2 - (    V2 - D)*HL1 - Q;

	H1.equals(0.0);

	//==========================================================================

	const T y7 = y_eq(x7);

	const T M6 = V7* y7 - (1 + V7 - D)* x6 - d;

	M6.equals(-0.5);

	//--------------------------------------------------------------------------

	const T HV7 = H_Vap(x7);

	const T HL6 = H_Liq(x6);

	const T H6 = V7*HV7 - (1 + V7 - D)*HL6 - Q;

	H6.equals(-0.0968047);

	//==========================================================================

	const T y3 = y_eq(x3);

	const T M2 = V3* y3 - (    V3 - D)* x2 - d;

	M2.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV3 = H_Vap(x3);

	const T HL2 = H_Liq(x2);

	const T H2 = V3*HV3 - (    V3 - D)*HL2 - Q;

	H2.equals(0.0);

	//==========================================================================

	const T y6 = y_eq(x6);

	const T M5 = V6* y6 - (1 + V6 - D)* x5 - d;

	M5.equals(-0.5);

	//--------------------------------------------------------------------------

	const T HV6 = H_Vap(x6);

	const T HL5 = H_Liq(x5);

	const T H5 = V6*HV6 - (1 + V6 - D)*HL5 - Q;

	H5.equals(-0.0968047);

	//==========================================================================

	const T y4 = y_eq(x4);

	const T M3 = V4* y4 - (    V4 - D)* x3 - d;

	M3.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV4 = H_Vap(x4);

	const T HL3 = H_Liq(x3);

	const T H3 = V4*HV4 - (    V4 - D)*HL3 - Q;

	H3.equals(0.0);

	//==========================================================================

	const T y5 = y_eq(x5);

	const T M4 = V5* y5 - (    V5 - D)* x4 - d;

	M4.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV5 = H_Vap(x5);

	const T HL4 = H_Liq(x4);

	const T H4 = V5*HV5 - (    V5 - D)*HL4 - Q;

	H4.equals(0.0);

	//cout << "V5: " << V5.compute_bounds() << endl;
	//cout << "x5: " << x5.compute_bounds() << endl;

	return;
}


//==============================================================================

void build(const problem<builder>* prob) {

	builder::reset();

	builder* box = prob->initial_box();

	prob->evaluate(box);

	delete[] box;

	ASSERT(builder::number_of_variables() == prob->number_of_variables())

	delete prob;

}

void dag_test(const problem<builder>* prob) {

	build(prob);

	// FIXME Recording occurence info should be part of the building procedure
	builder::record_occurence_info();

	expression_graph<interval> dag( builder::number_of_arguments(),
									builder::get_primitives(),
									builder::get_numeric_constants(),
									builder::get_rhs_of_constraints(),
									builder::get_initial_box());

	//builder::dbg_dump_type_of_primitives();
	builder::print_info(std::cout);
	builder::print_primitives(std::cout);
	builder::print_index_set(std::cout);
	builder::print_type1_common_subexpressions(std::cout);
	builder::print_type2_common_subexpressions(std::cout);
	builder::reset();

	for (int i=0; i<40; ++i) {
		dag.revise_all();
	}

	dag.show_variables(std::cout);

	std::cout << "Last value: " << dag.last_value() << std::endl;
}

void example_Hansen() {

	std::cout << "===============================================" << std::endl;
	std::cout << "Hansen\'s example" << std::endl;

	dag_test(new Hansen_example<builder> ());
}

void example_challenge() {

	std::cout << "===============================================" << std::endl;
	std::cout << "Neumaier\'s interval challange" << std::endl;

	dag_test(new Example_challange<builder> ());
}

void example_1() {

	std::cout << "===============================================" << std::endl;
	std::cout << "(x-1)/(x^2+2)" << std::endl;

	dag_test(new Example_1<builder> ());
}

void example_2() {

	std::cout << "===============================================" << std::endl;
	std::cout << "(x^2+x)/(16*x-9)" << std::endl;

	dag_test(new Example_2<builder> ());
}

void example_3() {

	std::cout << "===============================================" << std::endl;
	std::cout << "x^2+y^2=1" << std::endl;
	std::cout << "x^2-y  =0" << std::endl;
	dag_test(new Example_3<builder> ());
}

void example_Jacobsen() {

	std::cout << "===============================================" << std::endl;
	std::cout << "Jacobsen" << std::endl;

	dag_test(new Jacobsen<builder> ());
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

		std::cout << demangle(typeid(e).name()) << std::endl;
		std::cout << e.what() << std::endl;
	}

	example_Hansen();

	example_challenge();

	example_1();

	example_2();

	example_3();

	example_Jacobsen();

	return 0;
}
