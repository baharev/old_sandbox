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

#include <limits>
#include "affine_pair_iterator.hpp"
#include "diagnostics.hpp"

namespace asol {


affine_pair_iterator::affine_pair_iterator(const affine& x, const affine& y) :
	i    (x.noise_vars.begin()),
	i_end(x.noise_vars.end()),
	j    (y.noise_vars.begin()),
	j_end(y.noise_vars.end())
{
	x.dbg_consistency();
	y.dbg_consistency();

	k = i->index; // k = 0;

	x_k = i->coeff;
	y_k = j->coeff;
}

bool affine_pair_iterator::increment() {

	if (i == i_end && j == j_end) {

		return false;
	}

	const int int_max = std::numeric_limits<int>::max();

    const int i_index = (i != i_end) ? i->index : int_max;

    const int j_index = (j != j_end) ? j->index : int_max;

	if (i_index < j_index) {

		k = i_index;

		x_k = i->coeff;
		y_k = 0.0;

		++i;
	}
	else if (j_index < i_index ) {

		k = j_index;

		x_k = 0.0;
		y_k = j->coeff;

		++j;
	}
	else {

		ASSERT(i_index!=int_max && i_index==j_index);

		k = i_index;

		x_k = i->coeff;
		y_k = j->coeff;

		++i;
		++j;
	}

	return true;
}

}
