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

#include <ostream>
#include <sstream>
#include "printer.hpp"
#include "diagnostics.hpp"

using std::endl;

namespace asol {

printer::printer(std::ostream& os, const std::map<int,double>& num_const)
: out(os), numeric_const(num_const)
{

}

char printer::type(const int index) const {

	Map::const_iterator i = numeric_const.find(index);

	return (i==numeric_const.end()) ? 'v' : 'n' ;
}

const std::string printer::arg(const int index) const {

	std::ostringstream os;

	os << type(index) << index << std::flush;

	return os.str();
}

void printer::record_unary_primitive(int z, int x, const char* op) {

	out << 'v' << z << " = " << op << '(' << arg(x) << ')' << endl;
}

void printer::record_binary_primitive(int z, int x, int y, const char* op) {

	out << 'v' << z << " = " << arg(x) << ' ' << op << ' ' << arg(y) << endl;
}

void printer::addition(int z, int x, int y) {

	record_binary_primitive(z, x, y, "+");
}

void printer::substraction(int z, int x, int y) {

	record_binary_primitive(z, x, y, "-");
}

void printer::multiplication(int z, int x, int y) {

	record_binary_primitive(z, x, y, "*");
}

void printer::division(int z, int x, int y) {

	record_binary_primitive(z, x, y, "/");
}

void printer::square(int z, int x) {

	record_unary_primitive(z, x, "sqr");
}

void printer::exponential(int z, int x) {

	record_unary_primitive(z, x, "exp");
}

void printer::equality_constraint(int z, int x, double val) {

	out << "Eq: v" << z << " = rhs(" << x << ")" << endl << endl;
}

}
