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

#include <iosfwd>
#include <map>
#include <vector>
#include "typedefs.hpp"

namespace asol {

template<typename T> class primitive;

template <typename T>
class expression_graph {

public:

	typedef std::vector<primitive<T>*> PrimVector;
	typedef std::map<int,double> Map;

	expression_graph(int number_of_arguments,
			const PrimVector& primitives,
			const Map& numeric_constants,
			const BoundVector& initialbox);

	void evaluate_all();

	void revise_all();

	void show_variables(std::ostream& out) const;

	const T& last_value() const;

	void evaluate_primitive(int i);

	~expression_graph();

private:

	expression_graph(const expression_graph& );
	expression_graph& operator=(const expression_graph& );

	int v_size() const { return static_cast<int> (v.size()); }
	void set_variables();
	void set_non_variables(const int length);
	void set_numeric_consts(const int length);

	std::vector<T> v;
	const PrimVector primitives;
	const Map constants;
	const BoundVector initial_box;
};

}

#endif // EXPRESSION_GRAPH_HPP_
