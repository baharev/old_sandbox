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

#ifdef __GNUG__
#include <cxxabi.h>
#include <cstdlib>
#endif

#include "demangle.hpp"

using namespace std;

namespace asol {


#ifdef __GNUG__

const string name(const type_info& t) {

	const char* const name = t.name();

	int status = -4;

	char* res = abi::__cxa_demangle(name, NULL, NULL, &status);

	const char* const demangled_name = (status==0)?res:name;

	string ret_val(demangled_name);

	free(res);

	return ret_val;
}

#else

const string name(const type_info& t) {

	return string(t.name());
}

#endif

}
