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

#include <vector>
#include "interval.hpp"
#include "primitives.hpp"

using namespace std;

using namespace asol;

void Main() {

	vector<interval> v;

	v.push_back(interval(1, 2));
	v.push_back(interval(0, 1));
	v.push_back(interval(1, 3));
	v.push_back(interval(1, 9));

	primitive<interval>::set_vector(v);

	v.clear();

	primitive<interval>* p = new addition<interval>(2, 0, 1);

	primitive<interval>* q = new square<interval>(3, 2);

	p->evaluate();
	q->evaluate();

	p->revise();
	q->revise();

	delete p;
	delete q;

}

int main() {

	Main();

	return 0;
}
