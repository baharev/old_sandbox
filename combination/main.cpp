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

#include <iostream>
#include <iterator>
#include "box_generator.hpp"

using namespace std;
using namespace asol;

void Main() {

	box_generator::IVector vars;

	vars.push_back( interval( 0, 3) );
	vars.push_back( interval( 2, 5) );
	vars.push_back( interval( 1, 4) );

	box_generator::IntVector index_set;

	index_set.push_back(0);
	index_set.push_back(1);
	index_set.push_back(2);

	box_generator generator(vars, index_set, 3);

	while (generator.set_next()) {

		copy(vars.begin(), vars.end(), ostream_iterator<interval>(cout,"\t"));

		cout << endl;
	}
}

int main() {

	Main();

	return 0;
}
