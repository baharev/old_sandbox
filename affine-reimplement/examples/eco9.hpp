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

#ifndef ECO9_HPP_
#define ECO9_HPP_

#include "problem.hpp"

namespace asol {

template <typename T>
class eco9 : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	virtual int number_of_stored_solutions() const;

	virtual const DoubleArray2D solutions() const;

	enum { X1, X2, X3, X4, X5, X6, X7, X8, SIZE, SOLS = 16 };
};

}

#endif // ECO9_HPP_
