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
		} ,
		{
				0,
				1,
				0,

				352.7488454095784,

				0.569532026149012,
				1,
				0.3347082310792669,

				0,
				1,
				0,

				0.19268114050611462,
				0,
				0.9999879349507176,

				0.8064660238711381,
				0,
				1.315264521636581
		} ,
		{
				0,
				0,
				1,

				373.1568320866717,

				0.18184547527890899,
				2.053400617544793e-05,
				1,

				0,
				0,
				1,

				0.07313260941786259,
				0.6121091417029763,
				0,

				2.5790518917962277,
				12.01601195416843,
				0
		} ,
		{
				0.4866448919196974,
				0.5133551086632554,
				0,

				347.2449377945518,

				0.7744613312802209,
				0.90756658564835,
				0.5755221088713942,

				0.628365642368807,
				0.5656390580934878,
				0,

				-0.08656711994290278,
				0.08206297359757586,
				0.8901645350887561,

				0,
				0,
				0.43556401271443823
		} ,
		{
				0.89209110494315,
				0,
				0.1079088951027853,

				351.25709309000126,

				0.9110753265296915,
				0.7208575732934747,
				0.8727694642517312,

				0.9791628408390195,
				0,
				0.12363963168131802,

				-0.08516891861867104,
				0.4037650481402758,
				0.7040979767698728,

				0,
				0.6831984321505225,
				0
		} ,
		{
				0,
				0.6805730476636678,
				0.31942695282040123,

				346.772836726618,

				0.4368340308152552,
				0.6805763069950197,
				0.5368776097927617,

				0,
				0.9999952109244498,
				0.594971641569672,

				-0.31172743277827175,
				-0.19009526672799026,
				0.40501815479624514,

				0.3282190988788769,
				0,
				0
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
	x[U2] = T(-0.3207,    0.62   ); // FIXME What is going on???
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

	x1.less_than_or_equal_to(t1);

	x2.less_than_or_equal_to(t2);

	x3.less_than_or_equal_to(t3);

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
const DoubleArray2D Wilson16<T>::solutions() const {

	ASSERT2(n_vars==SIZE,"n_vars: "<<n_vars)

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	DoubleArray2D solution_vectors(SOLS);

	for (int i=0; i<SOLS; ++i) {

		const double* const x = sol[i];

		solution_vectors.at(i).assign(x, x + SIZE);
	}

	return solution_vectors;
}

template class Wilson16<builder>;

}
