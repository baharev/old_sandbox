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
#include "algorithm.hpp"
#include "envelope.hpp"
#include "problem.hpp"

using namespace std;
using namespace asol;

void example() {

	var::reset();

	const var x1(0.6, 0.7);
	const var x2(0.2, 0.3);
	const var x3(0.1, 0.2);

	const var y1(0.1, 0.4);
	const var y2(0.2, 0.4);
	const var y3(0.3, 0.5);

	var sx = x1 + x2 + x3;
	var sy = y1 + y2 + y3;

	sx.fix_at(1);
	sy.fix_at(1);

	const var z = (x2-1)*(y3-1)-(x3-x2)*y2;

	//const var z = x1*y1 + x1*y2 + x2*y2 + x3*y1;
	// [0.23, 0.76] -> [0.386667, 0.60]

	cout << "z: " << z << endl;

	var::dump_lp("debug.txt");
}

void example_Hansen() {

	var::reset();

	const var x(1, 10);
	const var y(1, 10);

	const var xy = x*y;

	const var numerator   = 5*x-4*sqr(y)+14*xy;
	const var denominator = sqr(x)+y+xy;

	cout << endl << "xy" << endl << xy << endl;
	cout << endl << "5x-4y^2+14xy" << endl << numerator << endl;
	cout << endl << "x^2+y+xy" << endl << denominator << endl;

	const var z = numerator/denominator;

	cout << endl << "z = (5x-4y^2+14xy)/(x^2+y+xy)" << endl << z << endl;

	// Base variables are not tightened, that is why z is not sharp

	//x.tighten_bounds();
	//y.tighten_bounds();

	//cout << endl << "x" << endl << x << endl;
	//cout << endl << "y" << endl << y << endl;
}

void example_challange() {

	var::reset();

	const var w(-0.9, -0.6);
	const var x(-0.1,  0.2);

	const var u = sqr(w) + sqr(x);

	const var y( 0.3,  0.7);
	const var z(-0.2,  0.1);

	const var v = sqr(y) + sqr(z);

	const var b(-1, 1);
	const var c(-1, 1);

	const var d = b*(x*y - w*z) + c*(x*z + w*y);

	const var a(7.0, 9.0);

	const var f = (a*(u-v)+2*d)/(u+v);

	cout << "f: " << f << endl;
}

void example_digression() {

	var::reset();

	enum { W, X, Y, Z, N, D, SIZE };

	double lb[SIZE]; double ub[SIZE];

	lb[W] = -0.9; ub[W] = -0.6;
	lb[X] = -0.1; ub[X] =  0.2;
	lb[Y] =  0.3; ub[Y] =  0.7;
	lb[Z] = -0.2; ub[Z] =  0.1;

	lb[N] = -100; ub[N] =  100;
	lb[D] = -100; ub[D] =  100;

	for (int i=1; i<=2; ++i) {

		var w(lb[W], ub[W]);
		var x(lb[X], ub[X]);
		var y(lb[Y], ub[Y]);
		var z(lb[Z], ub[Z]);

		var num = 2*(x*z+w*y);

		num.intersect(lb[N], ub[N]);

		var den = sqr(w)+sqr(x)+sqr(y)+sqr(z);

		den.intersect(lb[D], ub[D]);

		var s = num/den;

		cout << "s: " << s << endl;

		if (i==2)
			break;

		num.tighten_bounds();
		den.tighten_bounds();

		w.tighten_bounds();
		x.tighten_bounds();
		y.tighten_bounds();
		z.tighten_bounds();

		num.copy_bounds(lb[N], ub[N]);
		den.copy_bounds(lb[D], ub[D]);

		w.copy_bounds(lb[W], ub[W]);
		x.copy_bounds(lb[X], ub[X]);
		y.copy_bounds(lb[Y], ub[Y]);
		z.copy_bounds(lb[Z], ub[Z]);

		cout << "w: " << w << endl;
		cout << "x: " << x << endl;
		cout << "y: " << y << endl;
		cout << "z: " << z << endl;

		var::reset();
	}
}

void example_1() {

	var::reset();

	const var x(2.0, 4.0);

	const var y = (x-1)/(sqr(x)+2);

	cout << endl << "y = (x-1)/(x^2+2)" << endl << y << endl;

}

void example_2() {

	var::reset();

	const var x(1.0, 2.0);

	const var y = (sqr(x)+x)/(16*x-9);

	cout << endl << "y = (x^2+x)/(16*x-9)" << endl << y << endl;
}

void example_Wilson() {

	var::reset();

	const var x1(0, 1);
	const var x2(0, 1);
	const var x3(0, 1);

	var sx = x1+x2+x3;

	sx.fix_at(1);

	const var Lambda12(0.531930, 0.611472);
	const var Lambda13(0.169655, 0.183591);
	const var Lambda21(0.795033, 0.819308);
	const var Lambda23(6.17202e-06, 2.42304e-05);
	const var Lambda31(0.787092, 0.947967);
	const var Lambda32(0.277934, 0.406111);

    const var s1 =          x1 + Lambda12*x2 + Lambda13*x3;
    const var s2 = Lambda21*x1 +          x2 + Lambda23*x3;
    const var s3 = Lambda31*x1 + Lambda32*x2 +          x3;

    cout << "s1 : " << s1 << endl;
    cout << "s2 : " << s2 << endl;
    cout << "s3 : " << s3 << endl;

    const var t1 = x1/s1;
    const var t2 = x2/s2;
    const var t3 = x3/s3;

    cout << "t1 : " << t1 << endl;
    cout << "t2 : " << t2 << endl;
    cout << "t3 : " << t3 << endl;

    //t1.intersect(0, 1);
    //t2.intersect(0, 1);
    //t3.intersect(0, 1);

    const var u1 =          t1 + Lambda21*t2 + Lambda31*t3;
    const var u2 = Lambda12*t1 +          t2 + Lambda32*t3;
    const var u3 = Lambda13*t1 + Lambda23*t2 +          t3;

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

class Jacobsen : public problem {

	virtual int size() const;

	virtual interval* initial_box() const;

	virtual void evaluate(const var x[]) const;

	enum {
		X1, X2, X3, X4, X5, X6, X7, X8,
		v1, v2, v3, v4, v5, v6, v7, C, SIZE
	};
};

int Jacobsen::size() const {

	return SIZE;
}

interval* Jacobsen::initial_box() const {

	interval* v = new interval[SIZE];

	//      interval(0.93, 0.94);
	v[X1] = interval(1.0e-4, 1.0);

	for (int i=X2; i<=X8; ++i) {
		v[i] = interval(1.0e-4, 1.0);
	}

	for (int i=v1; i<=v7; ++i) {
		v[i] = interval(2.0, 4.0);
	}
	//     interval(0.50, 0.51);
	v[C] = interval(0.0, 1.12);

	return v;
}

void Jacobsen::evaluate(const var v[]) const {

	const var& x1 = v[X1];
	const var& x2 = v[X2];
	const var& x3 = v[X3];
	const var& x4 = v[X4];
	const var& x5 = v[X5];
	const var& x6 = v[X6];
	const var& x7 = v[X7];
	const var& x8 = v[X8];
	const var& V1 = v[v1];
	const var& V2 = v[v2];
	const var& V3 = v[v3];
	const var& V4 = v[v4];
	const var& V5 = v[v5];
	const var& V6 = v[v6];
	const var& V7 = v[v7];
	const var& D = v[C];

	const var y1 = y_eq(x1);

	var d = D*y1;

	var Lw = (V1 - D)*(0.6010 - 0.2806*y1);

	Lw.fix_at(0.96);

	//--------------------------------------------------------------------------

	const var HV1 = H_Vap(x1);

	const var HL0 = H_Liq(y1);

	var Q = V1*(HV1 - HL0) + D*HL0;

	//==========================================================================

	var M8 = (1.0-D)*x8 + d;

	M8.fix_at(0.5);

	//==========================================================================

	const var L7 = 4.0-D;

	var y8 = y_eq(x8);

	var M7 = 3.0*y8-L7*x7-d;

	M7.fix_at(-0.5);

	//--------------------------------------------------------------------------

	const var HV8 = H_Vap(x8);

	const var HL7 = H_Liq(x7);

	var H7 = 3*HV8 -L7*HL7 - Q;

	H7.fix_at(-0.0968047);

	//==========================================================================

	const var y2 = y_eq(x2);

	var M1 = V2* y2 - (    V2 - D)* x1 - d;

	M1.fix_at(0.0);

	//--------------------------------------------------------------------------

	const var HV2 = H_Vap(x2);

	const var HL1 = H_Liq(x1);

	var H1 = V2*HV2 - (    V2 - D)*HL1 - Q;

	H1.fix_at(0.0);

	//==========================================================================

	const var y7 = y_eq(x7);

	var M6 = V7* y7 - (1 + V7 - D)* x6 - d;

	M6.fix_at(-0.5);

	//--------------------------------------------------------------------------

	const var HV7 = H_Vap(x7);

	const var HL6 = H_Liq(x6);

	var H6 = V7*HV7 - (1 + V7 - D)*HL6 - Q;

	H6.fix_at(-0.0968047);

	//==========================================================================

	const var y3 = y_eq(x3);

	var M2 = V3* y3 - (    V3 - D)* x2 - d;

	M2.fix_at(0.0);

	//--------------------------------------------------------------------------

	const var HV3 = H_Vap(x3);

	const var HL2 = H_Liq(x2);

	var H2 = V3*HV3 - (    V3 - D)*HL2 - Q;

	H2.fix_at(0.0);

	//==========================================================================

	const var y6 = y_eq(x6);

	var M5 = V6* y6 - (1 + V6 - D)* x5 - d;

	M5.fix_at(-0.5);

	//--------------------------------------------------------------------------

	const var HV6 = H_Vap(x6);

	const var HL5 = H_Liq(x5);

	var H5 = V6*HV6 - (1 + V6 - D)*HL5 - Q;

	H5.fix_at(-0.0968047);

	//==========================================================================

	const var y4 = y_eq(x4);

	var M3 = V4* y4 - (    V4 - D)* x3 - d;

	M3.fix_at(0.0);

	//--------------------------------------------------------------------------

	const var HV4 = H_Vap(x4);

	const var HL3 = H_Liq(x3);

	var H3 = V4*HV4 - (    V4 - D)*HL3 - Q;

	H3.fix_at(0.0);

	//==========================================================================

	const var y5 = y_eq(x5);

	var M4 = V5* y5 - (    V5 - D)* x4 - d;

	M4.fix_at(0.0);

	//--------------------------------------------------------------------------

	const var HV5 = H_Vap(x5);

	const var HL4 = H_Liq(x4);

	var H4 = V5*HV5 - (    V5 - D)*HL4 - Q;

	H4.fix_at(0.0);

	//cout << "V5: " << V5.compute_bounds() << endl;
	//cout << "x5: " << x5.compute_bounds() << endl;

	return;
}

int Main() {

	algorithm a(new Jacobsen);

	a.run();

	return 0;

	example_Hansen();

	example_1();

	example_2();

	example_challange();

	example_digression();

	example_Wilson();

	return 0;
}

int main() {

	Main();

	return 0;
}
