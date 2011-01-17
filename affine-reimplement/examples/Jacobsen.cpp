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

#include "Jacobsen.hpp"
#include "builder.hpp"

namespace asol {

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

	y1.mark_as_common_subexpression();
	
	const T d = D*y1;

	d.mark_as_common_subexpression();

	const T Lw = (V1 - D)*(0.6010 - 0.2806*y1);

	Lw.equals(0.96);

	//--------------------------------------------------------------------------

	const T HV1 = H_Vap(x1);

	const T HL0 = H_Liq(y1);

	const T Q = V1*(HV1 - HL0) + D*HL0; // HL0*(D - V1) + HV1*V1 is worse...

	Q.mark_as_common_subexpression();

	//==========================================================================

	const T M8 = (1.0-D)*x8 + d;

	M8.equals(0.5);

	//==========================================================================

	const T L7 = 4.0-D;

	L7.mark_as_common_subexpression();

	const T y8 = y_eq(x8);

	const T M7 = 3.0*y8-L7*x7-d;

	M7.equals(-0.5);

	//--------------------------------------------------------------------------

	const T HV8 = H_Vap(x8);

	const T HL7 = H_Liq(x7);

	const T H7 = 3*HV8 -L7*HL7 - Q;

	H7.equals(-0.0968047);

	//==========================================================================

	const T L1 = V2 - D;

	L1.mark_as_common_subexpression();

	const T y2 = y_eq(x2);

	const T M1 = V2*y2 - L1*x1 - d;

	M1.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV2 = H_Vap(x2);

	const T HL1 = H_Liq(x1);

	const T H1 = V2*HV2 - L1*HL1 - Q;

	H1.equals(0.0);

	//==========================================================================

	const T L6 = 1 + V7 - D;

	L6.mark_as_common_subexpression();

	const T y7 = y_eq(x7);

	const T M6 = V7* y7 - L6*x6 - d;

	M6.equals(-0.5);

	//--------------------------------------------------------------------------

	const T HV7 = H_Vap(x7);

	const T HL6 = H_Liq(x6);

	const T H6 = V7*HV7 - L6*HL6 - Q;

	H6.equals(-0.0968047);

	//==========================================================================

	const T L2 = V3 - D;

	L2.mark_as_common_subexpression();

	const T y3 = y_eq(x3);

	const T M2 = V3*y3 - L2*x2 - d;

	M2.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV3 = H_Vap(x3);

	const T HL2 = H_Liq(x2);

	const T H2 = V3*HV3 - L2*HL2 - Q;

	H2.equals(0.0);

	//==========================================================================

	const T L5 = 1 + V6 - D;

	L5.mark_as_common_subexpression();

	const T y6 = y_eq(x6);

	const T M5 = V6*y6 - L5*x5 - d;

	M5.equals(-0.5);

	//--------------------------------------------------------------------------

	const T HV6 = H_Vap(x6);

	const T HL5 = H_Liq(x5);

	const T H5 = V6*HV6 - L5*HL5 - Q;

	H5.equals(-0.0968047);

	//==========================================================================

	const T L3 = V4 - D;

	L3.mark_as_common_subexpression();

	const T y4 = y_eq(x4);

	const T M3 = V4*y4 - L3*x3 - d;

	M3.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV4 = H_Vap(x4);

	const T HL3 = H_Liq(x3);

	const T H3 = V4*HV4 - L3*HL3 - Q;

	H3.equals(0.0);

	//==========================================================================

	const T L4 = V5 - D;

	L4.mark_as_common_subexpression();

	const T y5 = y_eq(x5);

	const T M4 = V5*y5 - L4*x4 - d;

	M4.equals(0.0);

	//--------------------------------------------------------------------------

	const T HV5 = H_Vap(x5);

	const T HL4 = H_Liq(x4);

	const T H4 = V5*HV5 - L4*HL4 - Q;

	H4.equals(0.0);

	return;
}

template class Jacobsen<builder>;

}
