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

#ifndef EXPRESSION_GRAPH_HPP_
#define EXPRESSION_GRAPH_HPP_

#include <vector>
#include "operations.hpp"
#include "primitives.hpp"

namespace asol {

template <typename T>
class expression_graph : private operations {

public:

	typedef std::vector<primitive*> Vector;

	expression_graph(int number_of_variables, const Vector& p){

		v.resize(number_of_variables);
		primitives = Vector(p);
	}

	~expression_graph() {

		clear_all();
	}

	void evaluate_primitive(int i) {

		primitives.at(i)->evaluate(this);
	}

private:

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	void clear_all() {

		v.clear();

		typedef Vector::iterator itr;

		for (itr i = primitives.begin(); i!=primitives.end(); ++i) {

			delete *i;
		}

		primitives.clear();
	}

	virtual void addition(int z, int x, int y) {

		v.at(z) = v.at(x) + v.at(y);
	}

	std::vector<T> v;

	std::vector<primitive*> primitives;
};

};

#endif // EXPRESSION_GRAPH_HPP_
