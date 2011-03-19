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

#include "eco9.hpp"
#include "builder.hpp"
#include "diagnostics.hpp"

namespace {


const int n_vars = 8;

const double sol[][n_vars] = {

		{
				-0.125,
				-0.125,
				-0.125,
				-0.125,
				-0.125,
				-0.125,
				-0.125,
				-0.125
		} ,
		{
				0.6455058742644363,
				0.6171284792730528,
				0.6728423906544725,
				-2.0910489241156136,
				0.14659704836498325,
				0.01246164440858569,
				-0.18730814164825882,
				-0.8161783711469645
		} ,
		{
				0.8528774451382914,
				0.12936397674410252,
				-1.1990131070180534,
				0.5776691420531274,
				-0.607113723146804,
				-0.27933102743596866,
				-0.010268933786665282,
				-0.4641837728495547
		} ,
		{
				-0.06341481490310709,
				-0.662094615365951,
				1.6233583590583476,
				0.12244932860834398,
				0.1833456972763919,
				-0.9509892630478123,
				-0.6064283617429337,
				-0.6462263298832794
		} ,
		{
				-0.039462627507042124,
				1.1335949725694425,
				-0.4726692536723607,
				-0.9487740868839563,
				-0.01972105278599304,
				0.14481050241442114,
				-0.3810670267752031,
				-0.41671142735930844
		} ,
		{
				0.5599685017999332,
				-0.7000568991841267,
				0.20027878701240168,
				-0.44304502788931566,
				0.016634108634041796,
				-0.311823976795595,
				-0.07712828296547564,
				-0.24482721061186385
		} ,
		{
				0.0221225575898636,
				0.6017681921211852,
				1.307916726919901,
				-1.2444837040428371,
				2.583056919137777,
				-0.27869129680155785,
				-1.8373702243416685,
				-2.1543191705826636
		} ,
		{
				1.5378459442063175,
				0.22412232503068727,
				-1.1781112699451606,
				-0.25855222000597883,
				0.3860831557558629,
				-1.4048614438733482,
				0.6026320549057931,
				-0.9091585460741732
		} ,
		{
				1.0,
				1.0,
				1.0,
				1.0,
				1.0,
				1.0,
				1.0,
				-8.0
		} ,
		{
				0.9384148149030948,
				1.4716040171553415,
				-0.28371746677282694,
				-0.18948365757618948,
				-2.512058521898249,
				1.0245553063204893,
				0.09813096738807803,
				-1.5474454595197378
		} ,
		{
				0.22949412888617277,
				-0.015268288045803186,
				-0.17961398483614852,
				2.5619999309891663,
				-0.8243815520764708,
				-0.7561196166930678,
				-0.7908882443801101,
				-1.22522237322567
		} ,
		{
				0.3150314981981559,
				1.2736491832961505,
				-0.0679422381075374,
				1.8096233127202164,
				-0.8180413710483373,
				2.859391721822124,
				-2.2871987980665462,
				-4.084513308814226
		} ,
		{
				0.9144626299033382,
				-0.3475078753536017,
				0.047325440232302216,
				2.2768132217947197,
				1.134284359245777,
				-2.720335706326677,
				0.09470013256151105,
				-2.3997422085494335
		} ,
		{
				1.6233833166991098,
				1.6249527705281592,
				0.6147114623961069,
				-1.0257825507946874,
				-2.044281748550139,
				-1.0303698082978936,
				2.268240601464076,
				-3.0308540434447315
		} ,
		{
				-0.6628459442061562,
				1.5452326219001988,
				-0.42465987634690994,
				0.28438628347294853,
				1.2958259866021251,
				-0.24651621655922518,
				-1.6915046895244799,
				-1.0999181653385015
		} ,
		{
				-0.7483833170716203,
				0.33996022059439557,
				0.6744903317165389,
				0.3384467001249095,
				-0.20281790314063947,
				-0.5361369262321265,
				-0.5356191005044688,
				-0.32994000548698893
		}
};

const int n_sol = sizeof(sol)/(n_vars*sizeof(double));

}

namespace asol {

template <typename T>
int eco9<T>::number_of_variables() const {

	ASSERT2(SIZE==n_vars,"SIZE, n_vars: "<<SIZE<<", "<<n_vars);
	return SIZE;
}

template <typename T>
T* eco9<T>::initial_box() const {

	T* x = new T[SIZE];

	for (int i=0; i<SIZE; ++i) {

		x[i] = T(-100.0, 100.0);
	}

	return x;
}

template <typename T>
void eco9<T>::evaluate(const T v[]) const {

	const T& x1 = v[X1];
	const T& x2 = v[X2];
	const T& x3 = v[X3];
	const T& x4 = v[X4];
	const T& x5 = v[X5];
	const T& x6 = v[X6];
	const T& x7 = v[X7];
	const T& x8 = v[X8];

	const T eq1 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;

	eq1.equals(-1);

	const T eq2 = x7 - x8*((7.0/8.0) - x1);

	eq2.equals(0.0);

	const T eq3 = x6 + x1*x7 - x8*((6.0/8.0) - x2);

	eq3.equals(0.0);

	const T eq4 = x5 + x1*x6 + x2*x7 - x8*((5.0/8.0) - x3);

	eq4.equals(0.0);

	const T eq5 = x4 + x1*x5 + x2*x6 + x3*x7 - x8*((4.0/8.0) - x4);

	eq5.equals(0.0);

	const T eq6 = x3*(1.0 + x6) + x4*(x1 + x7) + x2*x5 - x8*((3.0/8.0) - x5);

	eq6.equals(0.0);

	const T eq7 = x2 + x3*(x1 + x5) + x4*(x2 + x6) + x5*x7 - x8*((2.0/8.0) - x6);

	eq7.equals(0.0);

	const T eq8 = x1 + x2*(x1 + x3) + x4*(x3 + x5) + x6*(x5 + x7) - x8*((1.0/8.0) - x7);

	eq8.equals(0.0);

}

template <typename T>
int eco9<T>::number_of_stored_solutions() const {

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	return SOLS;
}

template <typename T>
const DoubleArray2D eco9<T>::solutions() const {

	ASSERT2(n_vars==SIZE,"n_vars: "<<n_vars)

	ASSERT2(SOLS==n_sol,"n_sol: "<<n_sol);

	DoubleArray2D solution_vectors(SOLS);

	for (int i=0; i<SOLS; ++i) {

		const double* const x = sol[i];

		solution_vectors.at(i).assign(x, x + SIZE);
	}

	return solution_vectors;
}

template class eco9<builder>;

}
