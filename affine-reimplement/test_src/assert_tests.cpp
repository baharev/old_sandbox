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
#include "diagnostics.hpp"
#include "builder.hpp"
#include "demangle.hpp"
#include "interval.hpp"

using namespace std;

namespace asol {

void assert_tests() {

	double lb(0), ub(2), x(1);

	ASSERT2(lb<=x && x<=ub, "lb, x, ub: "<<lb<<", "<<x<<", "<<ub);

	ASSERT(lb<=x && x<=ub);

	interval a(1, 3), b(-2, 5);

	interval c;

	//c.diameter();

	a/b;

	builder p;

	builder q;

	p+q;
}

void run_assert_test() {

	try {

		assert_tests();
	}
	catch (exception& e) {

		cout << demangle(typeid(e).name()) << endl;
		cout << e.what() << endl;
	}
}

}
