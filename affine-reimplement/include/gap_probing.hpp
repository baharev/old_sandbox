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

#ifndef GAP_PROBING_HPP_
#define GAP_PROBING_HPP_

#include <vector>

class interval;
template <typename T> class expression_graph;

namespace asol {

class gap_probing {

public:

	gap_probing(expression_graph<interval>& graph, interval* box, int length);

	// FIXME Hideous, pass the results in the graph instead?
	interval* contracted_box();

private:

	expression_graph<interval>& graph;

	const int n_vars;

	interval* box;

	std::vector<interval*> pending;

};

}

#endif // GAP_PROBING_HPP_
