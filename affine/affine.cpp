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
#include <limits>
#include <assert.h>
#include "affine.hpp"

using namespace std;

namespace asol {

int affine::max_used_index(0);

void affine::reset_max_used_index() {

	max_used_index = 0;
}

void affine::set_max_used_index(int idx) {

	max_used_index = idx;
}

int affine::get_max_used_index() {

	return max_used_index;
}

void affine::check_consistency() const {
	assert(range.inf() < range.sup());
	assert(n >= 2); // TODO Need for constant affine forms, n=1?
	assert(index[0] == 0);
	// TODO Is in LP? max used index?
}

void check_consistency(const affine& x, const affine& y) {
	x.check_consistency();
	y.check_consistency();
}

affine::affine(double lb, double ub) :
				index(new int[2]), value(new double[2]), n(2), range(interval(lb, ub))
{
	assert(lb < ub);

	index[0] = 0;
	index[1] = ++max_used_index;
	value[0] = range.midpoint();
	value[1] = range.radius();

	// TODO Push col into LP
}

affine::affine(const affine& x) :
				index(new int[x.n]), value(new double[x.n]), n(x.n), range(x.range)
{
	x.check_consistency();

	for (int i=0; i<n; ++i) {
		index[i] = x.index[i];
		value[i] = x.value[i];
	}
}

affine::affine(bool, int size) :
				index(new int[size]), value(new double[size]), n(0)
{

}

affine::~affine() {

	delete[] index;
	delete[] value;
}

std::ostream& operator<<(std::ostream& os, const affine& x) {

	os << "n: " << x.n << endl;

	for (int i=0; i<x.n; ++i) {
		os << x.index[i] << ": " << x.value[i] << endl;
	}

	os << "IA: " << x.range << endl;

	return os;
}

int affine::index_at(int i) const {

	assert(0<=i && i<=n);

	return i!=n ? index[i] : numeric_limits<int>::max();
}

bool next(const affine& x, const affine& y, int& idx, double& xi, double& yi) {

	static int i = 1;
	static int j = 1;

	if (idx==0) {
		i = j = 1;
	}

	if (i==x.n && j==y.n) {
		return false;
	}

	int x_ind = x.index_at(i);
	int y_ind = y.index_at(j);

	if (x_ind < y_ind) {
		assert(i<x.n);
		idx = x_ind;
		xi = x.value[i++];
		yi = 0.0;
	}
	else if (x_ind > y_ind) {
		assert(j<y.n);
		idx = y_ind;
		xi = 0.0;
		yi = y.value[j++];
	}
	else {
		assert(i<x.n && j<y.n);
		idx = x_ind;
		xi = x.value[i++];
		yi = y.value[j++];
	}

	return true;
}


const affine operator+(const affine& x, const affine& y) {

	check_consistency(x, y);

	const int size = x.n + y.n - 1;

	affine z(false, size);

	z.index[0] = 0;
	z.value[0] = x.value[0] + y.value[0];
	z.range = x.range + y.range;

	int k=1, index = 0;

	double x_i, y_i;

	while ( next(x, y, index, x_i, y_i) ) {

		assert (k<size);

		z.index[k] = index;
		z.value[k] = x_i + y_i;

		// TODO Compute mid and rad

		++k;
	}

	z.n = k;

	return z;
}



}
