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
#include "combination.hpp"
#include "interval.hpp"

using namespace std;
using namespace asol;

void Main() {

	std::vector<index_range> var;

	interval v[] = { interval(0,1), interval(2,3), interval(-3,1)};

	var.push_back(index_range(23,v  ));
	var.push_back(index_range( 4,v+1));
	var.push_back(index_range( 5,v+2));

	combination c(var, 3);

}

int main() {

	Main();

	return 0;
}