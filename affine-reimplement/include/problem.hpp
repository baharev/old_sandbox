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

#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

#include <vector>

namespace asol {

typedef std::vector<std::vector<double> > DoubleArray2D;

template <typename T>
class problem {

public:

	virtual int number_of_variables() const = 0;

	virtual T* initial_box() const = 0;

	virtual void evaluate(const T x[]) const = 0;

	virtual int number_of_stored_solutions() const { return 0; }

	virtual const DoubleArray2D solutions() const { return DoubleArray2D(); }

	virtual ~problem() { }

};

}

#endif // PROBLEM_HPP_
