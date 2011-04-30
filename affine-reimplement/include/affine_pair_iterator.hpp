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

#ifndef AFFINE_PAIR_ITERATOR_HPP_
#define AFFINE_PAIR_ITERATOR_HPP_

#include "affine.hpp"

namespace asol {

class affine_pair_iterator {

public:

	affine_pair_iterator(const affine& x, const affine& y);

	bool increment();

	int  index() const { return k; }
	double x_i() const { return x_k; }
	double y_i() const { return y_k; }

private:

	typedef std::vector<epsilon>::const_iterator itr;

	itr i;
	itr i_end;
	itr j;
	itr j_end;

	int k;
	double x_k;
	double y_k;
};

}

#endif // AFFINE_PAIR_ITERATOR_HPP_
