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

#include <algorithm>
#include <iostream>
#include <vector>
#include "splitting_strategy.hpp"
#include "diagnostics.hpp"
#include "interval.hpp"

using std::cout;
using std::endl;

namespace asol {

struct width {

	double operator()(const interval& x) { return x.diameter(); }
};

int find_max_diam_element(const interval* box, const int size) {

	double diameter[size];

	std::transform(box, box+size, diameter, width());

	int index = std::max_element(diameter, diameter+size) - diameter;

	return index;
}

int max_diam_selector::index_to_split(const interval* box) const {

	const int index = find_max_diam_element(box, n_vars);

	cout << "Splitting " << index << ", " << box[index] << endl;

	//ASSERT ( ! box[index].is_narrow(CONVERGENCE_TOL) ); // TODO Inconsistent with width

	return index;
}

Jacobsen_x1_D::Jacobsen_x1_D(const int n_vars) : splitting_strategy(n_vars) {

	ASSERT(n_vars == 16);
}

int Jacobsen_x1_D::index_to_split(const interval* box) const {

	double x1 = box[0].diameter();
	double D  = box[15].diameter();

	int    index = (x1 > D)? 0 : 15;
	double value = (x1 > D)? x1: D ;

	double V1 = box[8].diameter();

	if (V1 > value*2) {
		//index = 8;
	}

	return index;
}

eco9_sparsity::eco9_sparsity(const int n_vars) : splitting_strategy(n_vars) {

	ASSERT(n_vars == 8);
}

int eco9_sparsity::index_to_split(const interval* box) const {

	int index = widest_containing_zero(box);

	if (index == -1) {

		index = find_max_diam_element(box, n_vars);
	}

	cout << "Splitting " << index << ", " << box[index] << endl;

	return index;
}

int eco9_sparsity::widest_containing_zero(const interval* box) const {

	int index_set[] = { 0, 7, 1, 2, 3 };

	const int size = sizeof index_set / sizeof index_set[0];

	std::vector<int>      index;
	std::vector<interval> value;

	for (int i=0; i<size; ++i) {

		const interval& component =box[index_set[i]];

		if (component.contains(0) && !component.is_narrow()) {
			index.push_back(i);
			value.push_back(component);
		}
	}

	if (index.empty()) {

		return -1;
	}

	int selected = find_max_diam_element(&value.at(0), value.size());

	return index_set[index.at(selected)];
}

}
