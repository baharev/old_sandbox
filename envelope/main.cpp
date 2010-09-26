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

	var numerator   = 5*x-4*sqr(y)+14*xy;
	var denominator = sqr(x)+y+xy;

	cout << endl << "xy" << endl << xy << endl;
	cout << endl << "5x-4y^2+14xy" << endl << numerator << endl;
	cout << endl << "x^2+y+xy" << endl << denominator << endl;

	const var z = numerator/denominator;

	cout << endl << "z = (5x-4y^2+14xy)/(x^2+y+xy)" << endl << z << endl;
}


void example_1() {

	var x(2.0, 4.0);

	var y = (x-1)/(sqr(x)+2);

	cout << endl << "y = (x-1)/(x^2+2)" << endl << y << endl;

}

void example_2() {

	var x(1.0, 2.0);

	var y = (sqr(x)+x)/(16*x-9);

	cout << endl << "y = (x^2+x)/(16*x-9)" << endl << y << endl;
}

void example_Wilson() {

	var x1(0, 1);
	var x2(0, 1);
	var x3(0, 1);

	var sx = x1+x2+x3;

	sx.fix_at(1);

	var Lambda12(0.531930, 0.611472);
	var Lambda13(0.169655, 0.183591);
	var Lambda21(0.795033, 0.819308);
	var Lambda23(6.17202e-06, 2.42304e-05);
	var Lambda31(0.787092, 0.947967);
	var Lambda32(0.277934, 0.406111);

    var s1 =          x1 + Lambda12*x2 + Lambda13*x3;
    var s2 = Lambda21*x1 +          x2 + Lambda23*x3;
    var s3 = Lambda31*x1 + Lambda32*x2 +          x3;

    cout << "s1 : " << s1 << endl;
    cout << "s2 : " << s2 << endl;
    cout << "s3 : " << s3 << endl;

    var t1 = x1/s1;
    var t2 = x2/s2;
    var t3 = x3/s3;

    cout << "t1 : " << t1 << endl;
    cout << "t2 : " << t2 << endl;
    cout << "t3 : " << t3 << endl;

    //t1.intersect(0, 1);
    //t2.intersect(0, 1);
    //t3.intersect(0, 1);

    var u1 =          t1 + Lambda21*t2 + Lambda31*t3;
    var u2 = Lambda12*t1 +          t2 + Lambda32*t3;
    var u3 = Lambda13*t1 + Lambda23*t2 +          t3;

    cout << "u1 : " << u1 << endl;
    cout << "u2 : " << u2 << endl;
    cout << "u3 : " << u3 << endl;
/*
    var one(1, 1);

    var rs1 = one/s1;
    var rs2 = one/s2;
    var rs3 = one/s3;

    cout << "rs1 : " << rs1 << endl;
    cout << "rs2 : " << rs2 << endl;
    cout << "rs3 : " << rs3 << endl;

    var ru1 =           rs1 *x1 + (Lambda21*rs2)*x2 + (Lambda31*rs3)*x3;
    var ru2 = (Lambda12*rs1)*x1 +           rs2 *x2 + (Lambda32*rs3)*x3;
    var ru3 = (Lambda13*rs1)*x1 + (Lambda23*rs2)*x2 +           rs3 *x3;

    cout << "u1 : " << ru1 << endl;
    cout << "u2 : " << ru2 << endl;
    cout << "u3 : " << ru3 << endl;
*/
}

int main() {

	example_Wilson();

	return 0;
}
