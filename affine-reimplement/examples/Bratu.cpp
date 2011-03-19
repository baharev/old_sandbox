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

#include "Bratu.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"

namespace {

const int n_vars = 30;

const double sol[][n_vars] = {

		{
				0.017199401154330288,
				0.033340167381367754,
				0.04840507286011018,
				0.062377787155251226,
				0.07524294228122237,
				0.08698619727310501,
				0.09759429971867178,
				0.10705514371524098,
				0.11535782373205623,
				0.1224926838835868,
				0.12845136215136943,
				0.13322682913149778,
				0.13681342093116028,
				0.1392068658901154,
				0.14040430486090616,
				0.14040430484405003,
				0.1392068658403695,
				0.1368134208509439,
				0.13322682902468552,
				0.12845136202306567,
				0.1224926837398296,
				0.11535782357947057,
				0.10705514356066359,
				0.09759429956877351,
				0.0869861971340378,
				0.07524294215831821,
				0.062377787052799984,
				0.048405072781240635,
				0.033340167328053845,
				0.017199401127540377
		} ,

		{
				0.34882307255410716,
				0.6961712248205117,
				1.0414319085648864,
				1.3837443393370759,
				1.7219050397022702,
				2.05424349514289,
				2.3784644331821116,
				2.6914592276376292,
				2.9891021089473138,
				3.2660708639315086,
				3.5157678588365058,
				3.730457829555925,
				3.9017573792480538,
				4.02155917196816,
				4.083308913318545,
				4.08330891331809,
				4.021559171966856,
				3.9017573792460585,
				3.730457829553451,
				3.5157678588337706,
				3.266070863928702,
				2.9891021089445817,
				2.6914592276350753,
				2.378464433179805,
				2.0542434951408732,
				1.721905039700569,
				1.3837443393357058,
				1.041431908563855,
				0.6961712248198229,
				0.34882307255376244
		}

};

const int n_sol = sizeof(sol)/(n_vars*sizeof(double));

}

namespace asol {

template <typename T>
int Bratu<T>::number_of_variables() const {

	return SIZE;
}

template <typename T>
T* Bratu<T>::initial_box() const {

	T* x = new T[SIZE];

	for (int i=0; i<SIZE; ++i) {

		x[i] = T(-10.0, 10.0);
	}

	return x;
}

template <typename T>
void Bratu<T>::evaluate(const T x[]) const {

	const double c1 = -2.0, c2 = 1.04058273e-03;

	T eq_0 = x[1] + c1*x[0] + c2*exp(x[0]);

	eq_0.equals(0.0);

	for (int i=1; i<SIZE-1; ++i) {

		T eq_i = x[i+1] + c1*x[i] + c2*exp(x[i]) + x[i-1];

		eq_i.equals(0.0);
	}

	const int N = SIZE-1;

	T eq_N = x[N-1] + c1*x[N] + c2*exp(x[N]);

	eq_N.equals(0.0);
}

template <typename T>
int Bratu<T>::number_of_stored_solutions() const {

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	return SOLS;
}

template <typename T>
const DoubleArray2D Bratu<T>::solutions() const {

	ASSERT2(n_vars==SIZE,"n_vars: "<<n_vars)

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	DoubleArray2D solution_vectors(SOLS);

	for (int i=0; i<SOLS; ++i) {

		const double* const x = sol[i];

		solution_vectors.at(i).assign(x, x + SIZE);
	}

	return solution_vectors;
}

template class Bratu<builder>;

}
