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
#include "Bratu.hpp"
#include "eco9.hpp"
#include "Wilson16.hpp"
#include "expression_graph_test.hpp"
#include "builder.hpp"
#include "interval.hpp"
#include "search_procedure.hpp"

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

void example_Bratu() {

	cout << "###############################################" << endl;
	cout << "Bratu " << endl;

	test_probing_on_initial_box(new Bratu<builder> ());
}

void Bratu_solutions() {

	cout << "###############################################" << endl;
	cout << "Bratu solutions revise" << endl;

	test_solutions_revise(new Bratu<builder> ());
}

void Bratu_solutions_revise2() {

	cout << "###############################################" << endl;
	cout << "Bratu solutions revise2" << endl;

	test_solutions_revise2(new Bratu<builder> ());
}

void Jacobsen_solutions() {

	cout << "###############################################" << endl;
	cout << "Jacobsen solutions revise_all" << endl;

	test_solutions_revise(new Jacobsen<builder> ());
}

void Jacobsen_solutions_iterative_revise() {

	cout << "###############################################" << endl;
	cout << "Jacobsen solutions with iterative revision" << endl;

	test_solutions_iterative_revise(new Jacobsen<builder> ());
}

void Jacobsen_solutions_probing() {

	cout << "###############################################" << endl;
	cout << "Jacobsen solutions probing" << endl;

	test_solutions_probing(new Jacobsen<builder> ());
}

void eco9_solutions_iterative_revise() {

	cout << "###############################################" << endl;
	cout << "eco9 solutions iterative revise" << endl;

	test_solutions_iterative_revise(new eco9<builder> ());
}

void Wilson16_solutions_probing() {

	cout << "###############################################" << endl;
	cout << "Wilson16 solutions probing" << endl;

	dag_test(new Wilson16<builder> ());

	// FIXME Reconsider this index set stuff
	//extern void print_index_set(const problem<builder>* );
	//print_index_set(new Wilson16<builder> ());

	test_solutions_probing(new Wilson16<builder> ());
}

void eco9_extended_division_test() {

	cout << "###############################################" << endl;
	cout << "Extended division test on eco9" << endl;

	interval x[8];

	for (int i=0; i<8; ++i) {
		x[i] = interval(-100, 100);
	}

	x[6] = interval( 0.5, 100);
	x[7] = interval(-100,  -7);

	double sol[8] = { 1, 1, 1, 1, 1, 1, 1, -8 };

	extended_division_test(new eco9<builder> (), x, sol, 8);
}

void eco9_gap_probing() {

	cout << "###############################################" << endl;
	cout << "Gap probing on eco9" << endl;

	interval x[8];

	for (int i=0; i<8; ++i) {
		x[i] = interval(-100, 100);
	}

	x[6] = interval( 0.5, 100);
	x[7] = interval(-100,  -7);

	double sol[8] = { 1, 1, 1, 1, 1, 1, 1, -8 };

	gap_probing_test(new eco9<builder> (), x, sol, 8);
}

void index_recorder_test() {

	cout << "###############################################" << endl;
	cout << "Jacobsen index recorder test" << endl;

	index_recorder_test(new Jacobsen<builder> ());
}

void run_examples() {

	example_Hansen();

	example_challenge();

	example_1();

	example_2();

	example_3();

	example_Jacobsen();

	example_Bratu();

	Bratu_solutions();

	Bratu_solutions_revise2();

	Jacobsen_solutions();

	Jacobsen_solutions_iterative_revise();

	Jacobsen_solutions_probing();

	eco9_solutions_iterative_revise();

	Wilson16_solutions_probing();

	eco9_extended_division_test();

	eco9_gap_probing();

	builder::release();
}

void run_search_procedure() {

	search_procedure algorithm(new eco9<builder> ());

	algorithm.run();
}

}
