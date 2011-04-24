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
#include <functional>
#include <vector>
#include "affine.hpp"
#include "builder.hpp"
#include "interval.hpp"
#include "primitives.hpp"
#include "recorder.hpp"

using namespace std;

namespace {

typedef asol::primitive<asol::builder> Primitive;

typedef std::vector<Primitive*> Vector;

}

namespace asol {

template <typename T>
class converter : public recorder {

public:

	converter(const Vector& vec) {

		for_each(vec.begin(), vec.end(), bind2nd(mem_fun(&Primitive::record), this));
	}

	const vector<primitive<T>*>& result() const {

		return v;
	}

	~converter();

private:

	converter(const converter& );
	converter& operator=(const converter& );

	virtual void addition(int z, int x, int y) {

		v.push_back(new asol::addition<T>(z, x, y));
	}

	virtual void substraction(int z, int x, int y) {

		v.push_back(new asol::substraction<T>(z, x, y));
	}

	virtual void multiplication(int z, int x, int y) {

		v.push_back(new asol::multiplication<T>(z, x, y));
	}

	virtual void division(int z, int x, int y) {

		v.push_back(new asol::division<T>(z, x, y));
	}

	virtual void square(int z, int x) {

		v.push_back(new asol::square<T>(z, x));
	}

	virtual void exponential(int z, int x) {

		v.push_back(new asol::exponential<T>(z, x));
	}

	virtual void logarithm(int z, int x) {

		v.push_back(new asol::logarithm<T>(z, x));
	}

	virtual void equality_constraint(int z, int x, double val) {

		v.push_back(new asol::equality_constraint<T>(z, x, val));
	}

	virtual void common_subexpression(int z, int x) {

		v.push_back(new asol::common_subexpression<T>(z, x));
	}

	virtual void less_than_or_equal_to(int z, int x) {

		v.push_back(new asol::less_than_or_equal_to<T>(z, x));
	}

	vector<primitive<T>*> v;
};

template <typename T>
converter<T>::~converter() {
	// Dtor is out-of-line to make the compiler shut-up
}

template <typename T>
const vector<primitive<T>*> convert(const Vector& vec) {

	converter<T> conv(vec);

	return conv.result();
}

template
const vector<primitive<interval>*> convert(const Vector& );

template
const vector<primitive<affine>*> convert(const Vector& );

}
