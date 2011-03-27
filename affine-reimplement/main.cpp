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

#include <string>
#include "assert_tests.hpp"
#include "box_generator_tests.hpp"
#include "diagnostics.hpp"
#include "examples.hpp"

using std::string;
using namespace asol;


namespace {

const string SIMPLE_TESTS = "simple_tests";
const string SEARCH_PROC  = "search_procedure";

}

void simple_tests() {

	run_assert_test();

	run_box_generator_test();

	show_Jacobsen_sparsity();

	run_examples();
}

void search_procedure() {

	run_search_procedure(); // TODO Leaks, builder::release() is NOT called

}

int main(int argc, const char* argv[]) {

	ASSERT2(argc==2,"provide command line arguments");

	if (argv[1]==SIMPLE_TESTS) {

		simple_tests();
	}
	else if (argv[1]==SEARCH_PROC) {

		search_procedure();
	}
	else {

		ASSERT2(false,"command line argument not recognized: "<<argv[1]);
	}

	return 0;
}
