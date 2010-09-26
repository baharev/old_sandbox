//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2010 Ali Baharev
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
#include "envelope.hpp"

using namespace std;
using namespace asol;

void dummy_1() {

	var x(1, 2);

	var y(2, 4);

	var u = x+y;

	u.fix_at(6);

	cout << "u: " << u << endl;

	var v = x*y;

	v.fix_at(8);

	var w = u + v;

	var::dump_lp("debug.txt");
}

void example() {

	var x1(0.6, 0.7);
	var x2(0.2, 0.3);
	var x3(0.1, 0.2);

	var y1(0.1, 0.4);
	var y2(0.2, 0.4);
	var y3(0.3, 0.5);

	var sx = x1 + x2 + x3;
	var sy = y1 + y2 + y3;

	sx.fix_at(1);
	sy.fix_at(1);

	var z = (x2-1)*(y3-1)-(x3-x2)*y2;

	//var z = x1*y1 + x1*y2 + x2*y2 + x3*y1;
	// [0.23, 0.76] -> [0.386667, 0.60]

	cout << "z: " << z << endl;

	var::dump_lp("debug.txt");
}

void example_Hansen() {

	const var x(1.0, 10.0);
	const var y(1.0, 10.0);

	const var xy = x*y;

	const var numerator   = 5*x-4*sqr(y)+14*xy;
	const var denominator = sqr(x)+y+xy;

	cout << endl << "xy" << endl << xy << endl;
	cout << endl << "5x-4y^2+14xy" << endl << numerator << endl;
	cout << endl << "x^2+y+xy" << endl << denominator << endl;

	const var z = numerator/denominator;

	cout << endl << "z = (5x-4y^2+14xy)/(x^2+y+xy)" << endl << z << endl;
}

int main() {

	example_Hansen();

	return 0;
}
