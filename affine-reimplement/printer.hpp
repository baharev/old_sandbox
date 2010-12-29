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
#include <string>
#include "recorder.hpp"
#include "typedefs.hpp"

namespace asol {

class printer : public recorder {

public:

	printer(std::ostream& os, const PairVector& numeric_const);

private:

	virtual void record(const addition*            );
	virtual void record(const substraction*        );
	virtual void record(const multiplication*      );
	virtual void record(const division*            );
	virtual void record(const square*              );
	virtual void record(const exponential*         );
	virtual void record(const equality_constraint* );

	void record_unary_primitive(const unary_primitive* p, const char* op);
	void record_binary_primitive(const binary_primitive* p, const char* op);

	void dbg_check_if_sorted();
	int numeric_const_size() const;
	const std::string arg(const int index) const;
	char type(const int index) const;

	std::ostream& out;
	const PairVector& numeric_const;
};

}

#endif // PRINTER_HPP_
