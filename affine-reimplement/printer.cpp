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
#include "primitives.hpp"
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

void printer::record_unary_primitive(const unary_primitive* p, const char* op) {

	out << 'v' << p->z << " = " << op << '(' << arg(p->x) << ')' << endl;
}

void printer::record_binary_primitive(const binary_primitive* p, const char* op) {

	out << 'v' << p->z << " = " << arg(p->x) << ' ' << op << ' ' << arg(p->y) << endl;
}

void printer::record(const addition* p) {

	record_binary_primitive(p, "+");
}

void printer::record(const substraction* p) {

	record_binary_primitive(p, "-");
}

void printer::record(const multiplication* p) {

	record_binary_primitive(p, "*");
}

void printer::record(const division* p) {

	record_binary_primitive(p, "/");
}

void printer::record(const square* p) {

	record_unary_primitive(p, "sqr");
}

void printer::record(const exponential* p) {

	record_unary_primitive(p, "exp");
}

void printer::record(const equality_constraint* p) {

	out << "Eq: v" << p->z << " = rhs(" << p->x << ")" << endl << endl;
}

}
