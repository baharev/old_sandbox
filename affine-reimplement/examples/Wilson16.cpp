//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
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

#include "Wilson16.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"

namespace {

const int n_vars = 16;

const double sol[][n_vars] = {

		{
				0.232307,
				0.540553,
				0.22714,
				345.798,
				0.573705,
				0.728909,
				0.593567,
				0.404924,
				0.741592,
				0.382669,
				-0.327491,
				-0.0889876,
				 0.546717,
				0.0,
				0.0,
				0.0
		} ,
		{
				1,
				0,
				0,
				351.457,
				1,
				0.8079530807319288,
				0.8580310941724832,
				1,
				0,
				0,
				0,
				0.43253649084993007,
				0.8240150141906243,
				0,
				0.6043896455041708,
				0.14520011013979883
		}
};

const int n_sol = sizeof(sol)/(n_vars*sizeof(double));

}

namespace asol {

template <typename T>
int Wilson16<T>::number_of_variables() const {

	ASSERT2(SIZE==n_vars,"SIZE, n_vars: "<<SIZE<<", "<<n_vars);
	return SIZE;
}

template <typename T>
T* Wilson16<T>::initial_box() const {

	T* x = new T[SIZE];

	x[X1] = T(0, 1);
	x[X2] = T(0, 1);
	x[X3] = T(0, 1);

	x[t]  = T(330, 380);

	x[S1] = T(0.1750   , 1);
	x[S2] = T(1.828e-05, 1);
	x[S3] = T(0.3347   , 1);

	x[T1] = T(0, 1);
	x[T2] = T(0, 1);
	x[T3] = T(0, 1);

	x[U1] = T(-0.6536,    0.1927);
	x[U2] = T(-0.3207,    0.4326);
	x[U3] = T(-0.0003847, 0.99999);

	x[LN_K1] = T(-0.1353, 2.581);
	x[LN_K2] = T(-0.1416, 12.02);
	x[LN_K3] = T(-0.9722, 1.316);

	return x;
}

template <typename U>
void Wilson16<U>::evaluate(const U v[]) const {

	const U& x1 = v[X1];
	const U& x2 = v[X2];
	const U& x3 = v[X3];

	const U& T  = v[t];

	const U& s1 = v[S1];
	const U& s2 = v[S2];
	const U& s3 = v[S3];

	const U& t1 = v[T1];
	const U& t2 = v[T2];
	const U& t3 = v[T3];

	const U& u1 = v[U1];
	const U& u2 = v[U2];
	const U& u3 = v[U3];

	const U& ln_K1 = v[LN_K1];
	const U& ln_K2 = v[LN_K2];
	const U& ln_K3 = v[LN_K3];

	const U rT = (1.0/1.9858775)/T;

	rT.mark_as_common_subexpression();

	const U Lambda12 = (89.57/58.39)*exp(-694.0825*rT);
	Lambda12.mark_as_common_subexpression();

	const U Lambda13 = (18.05/58.39)*exp(-393.1971*rT);
	Lambda13.mark_as_common_subexpression();

	const U Lambda21 = (58.39/89.57)*exp(149.7978*rT);
	Lambda21.mark_as_common_subexpression();

	const U Lambda23 = (18.05/89.57)*exp(-6811.3433*rT);
	Lambda23.mark_as_common_subexpression();

	const U Lambda31 = (58.39/18.05)*exp(-926.263*rT);
	Lambda31.mark_as_common_subexpression();

	const U Lambda32 = (89.57/18.05)*exp(-1888.8509*rT);
	Lambda32.mark_as_common_subexpression();

	const U s1_con = x1          + x2*Lambda12 + x3*Lambda13 - s1;
    s1_con.equals(0);

    const U s2_con = x1*Lambda21 + x2          + x3*Lambda23 - s2;
    s2_con.equals(0);

    const U s3_con = x1*Lambda31 + x2*Lambda32 + x3          - s3;
    s3_con.equals(0);

    const U t1_con = x1/s1 - t1;
    t1_con.equals(0);

    const U t2_con = x2/s2 - t2;
    t2_con.equals(0);

    const U t3_con = x3/s3 - t3;
    t3_con.equals(0);

    const U u1_con =          t1 + Lambda21*t2 + Lambda31*t3 + u1;
    u1_con.equals(1);

    const U u2_con = Lambda12*t1 +          t2 + Lambda32*t3 + u2;
    u2_con.equals(1);

    const U u3_con = Lambda13*t1 + Lambda23*t2 +          t3 + u3;
    u3_con.equals(1);

    const U ln_K1_con = 3667.704902/(-46.976 + T) + log(s1) + ln_K1 - u1;
    ln_K1_con.equals(12.0457125667196);

    const U ln_K2_con = 2904.342681/(-51.191 + T) + log(s2) + ln_K2 - u2;
    ln_K2_con.equals(9.63112956671963);

    const U ln_K3_con = 3984.922839/(-39.734 + T) + log(s3) + ln_K3 - u3;
    ln_K3_con.equals(11.9515595667196);

    const U eq_1 = x1*ln_K1;
    eq_1.equals(0);

    const U eq_2 = x2*ln_K2;
    eq_2.equals(0);

    const U eq_3 = x3*ln_K3;
    eq_3.equals(0);

    const U sum_x = x1 + x2 + x3;
    sum_x.equals(1);
}

template <typename T>
int Wilson16<T>::number_of_stored_solutions() const {

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	return SOLS;
}

template <typename T>
const double* Wilson16<T>::solution(int i) const {

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	ASSERT2(0<=i&&i<SOLS,"i, SOLS: "<< i<<", "<<SOLS);

	return sol[i];
}

template class Wilson16<builder>;

}
