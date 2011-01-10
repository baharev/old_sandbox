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
using namespace std;

class Wheel { };

class Car {

public:

	void process(const vector<Wheel>& wheel) {

		for_each(wheel.begin(), wheel.end(), bind1st(mem_fun(&Car::put), this));
		for_each(wheel.begin(), wheel.end(), bind1st(mem_fun1_t<void,Car,Wheel>(&Car::put), this));
	}

private:

	void put(const Wheel w) { }
};

int main() {

	vector<Wheel> w(4);

	Car c;

	c.process(w);

	return 0;
}
