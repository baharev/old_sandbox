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

#ifndef DIAGNOSTICS_HPP_
#define DIAGNOSTICS_HPP_

#include <stdexcept>
#include <sstream>

// TODO What is a macro in a namespace anyway?

#ifdef ASOL_ENABLE_ASSERTS
#ifdef __GNUG__
#define FUNCTION_ __PRETTY_FUNCTION__
#else
#define FUNCTION_ __FUNCTION__
#endif

// TODO How can I eliminate code duplication in a macro without another macro?

#define ASSERT(condition) { \
	if (!(condition)) { \
		std::ostringstream os__; \
		os__ << "Assertion failed: " << #condition << ", "; \
		os__ << FUNCTION_ << " at "<< __FILE__ << ':' << __LINE__ << std::flush; \
		throw std::logic_error(os__.str()); \
	} \
}

#define ASSERT2(condition, message) { \
	if (!(condition)) { \
		std::ostringstream os__; \
		os__ << "Assertion failed: " << #condition << "; "; \
		os__ << message << "; " ; \
		os__ << FUNCTION_ << " at "<< __FILE__ << ':' << __LINE__ << std::flush; \
		throw std::logic_error(os__.str()); \
	} \
}
#elif defined ASOL_DISABLE_ASSERTS

#define ASSERT(condition) (void) 0;

#define ASSERT2(condition, message) (void) 0;

#else

#error ASSERT should be enabled or disabled

#endif

#endif // DIAGNOSTICS_HPP_
