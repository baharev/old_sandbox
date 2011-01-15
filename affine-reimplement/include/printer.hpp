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

#ifndef PRINTER_HPP_
#define PRINTER_HPP_

#include <iosfwd>
#include <map>
#include <string>
#include "recorder.hpp"

namespace asol {

class printer : public recorder {

public:

	printer(std::ostream& os, const std::map<int,double>& numeric_const);

private:

	virtual void addition      (int z, int x, int y);
	virtual void substraction  (int z, int x, int y);
	virtual void multiplication(int z, int x, int y);
	virtual void division      (int z, int x, int y);

	virtual void square     (int z, int x);
	virtual void exponential(int z, int x);

	virtual void equality_constraint(int z, int x, double val);

	void record_unary_primitive(int z, int x, const char* op);
	void record_binary_primitive(int z, int x, int y, const char* op);

	const std::string arg(const int index) const;
	char type(const int index) const;

	typedef std::map<int,double> Map;

	std::ostream& out;
	const Map& numeric_const;
};

}

#endif // PRINTER_HPP_
