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

#ifndef EXAMPLE_1_HPP_
#define EXAMPLE_1_HPP_

#include "problem.hpp"

namespace asol {

template <typename T>
class Example_1 : public problem<T> {

private:

	virtual int number_of_variables() const;

	virtual T* initial_box() const;

	virtual void evaluate(const T x[]) const;

	enum { X, SIZE };

};

}

#endif // EXAMPLE_1_HPP_
