//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
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

#include "lp_impl.hpp"

namespace asol {

lp_impl::lp_impl() {

	lp = glp_create_prob();

	parm = new glp_smcp;

	glp_init_smcp(parm);

	init();
}

lp_impl::~lp_impl() {

	delete parm;

	glp_delete_prob(lp);
}

void lp_impl::free_environment() {

	glp_free_env();
}

void lp_impl::init() {

	glp_set_obj_dir(lp, GLP_MIN);

	parm->presolve = GLP_OFF;

	parm->msg_lev = GLP_MSG_ON;

	//parm->meth = GLP_DUAL;
}

}
