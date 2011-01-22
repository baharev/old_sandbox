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
#include "Challenge.hpp"
#include "Example_1.hpp"
#include "Example_2.hpp"
#include "Example_3.hpp"
#include "Hansen.hpp"
#include "Jacobsen.hpp"
#include "expression_graph_test.hpp"
#include "builder.hpp"

using namespace std;

namespace asol {

void example_Hansen() {

	cout << "###############################################" << endl;
	cout << "Hansen\'s example" << endl;

	dag_test(new Hansen_example<builder> ());
}

void example_challenge() {

	cout << "###############################################" << endl;
	cout << "Neumaier\'s interval challange" << endl;

	dag_test(new Example_challange<builder> ());
}

void example_1() {

	cout << "###############################################" << endl;
	cout << "(x-1)/(x^2+2)" << endl;

	dag_test(new Example_1<builder> ());
}

void example_2() {

	cout << "###############################################" << endl;
	cout << "(x^2+x)/(16*x-9)" << endl;

	dag_test(new Example_2<builder> ());
}

void example_3() {

	cout << "###############################################" << endl;
	cout << "x^2+y^2=1" << endl;
	cout << "x^2-y  =0" << endl;
	dag_test(new Example_3<builder> ());
}

void example_Jacobsen() {

	cout << "###############################################" << endl;
	cout << "Jacobsen" << endl;

	dag_test(new Jacobsen<builder> ());
}

void show_Jacobsen_sparsity() {

	cout << "###############################################" << endl;
	cout << "Jacobsen sparsity" << endl;

	print_sparsity(new Jacobsen<builder> ());
}

void test_system_of_equations() {

	cout << "###############################################" << endl;
	cout << "Jacobsen system of equations" << endl;

	test_system_of_equations(new Jacobsen<builder> ());
}

void test_directed_revision() {

	cout << "###############################################" << endl;
	cout << "Jacobsen directed revision" << endl;

	test_directed_revision(new Jacobsen<builder> ());
}

void test_probing_Jacobsen() {

	cout << "###############################################" << endl;
	cout << "Jacobsen probing" << endl;

	test_probing_Jacobsen(new Jacobsen<builder> ());
}

void run_examples() {

	example_Hansen();

	example_challenge();

	example_1();

	example_2();

	example_3();

	example_Jacobsen();

	test_system_of_equations();

	test_directed_revision();

	test_probing_Jacobsen();

	builder::release();
}

}
