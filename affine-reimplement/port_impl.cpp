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

#include <iostream>
#include "port_impl.hpp"
#include "diagnostics.hpp"

extern "C" {

// Defined in LINPAS.f
void linpas_(double* A, int* M, int* N, int* IA, double* B, double* C, double* X, int* MAXITR, double* CTX, int* , double* SIMP, int* ISIMP, int* IE, int* ERRCOD);

// Defined at the bottom of the file
void cprint_(double* A, int* M, int* N, int* IA, double* B, double* C, double* X, double* CTX, int* IS, double* SIMP, int* ISIMP, int* IE, int* ITER, int* IPTG, int* IAG, int* IAS, double* U, int* DONE);

}

namespace {

int itr_tmp; // callback function cprint_ updates this at each iteration
}

namespace asol {

port_impl::port_impl() {

	A  = 0;
	M  = 0;
	N  = 0;
	IA = 0;
	B  = 0;
	C  = 0;
	X  = 0;
	MAXITR = 0;
	CTX    = 0;
	IS     = 0;
	SIMP   = 0;
	ISIMP  = 0;
	IE     = 0;
	IERR   = 0;

	rows_added   = 0;
	slacks_added = 0;

	previous_itr_count = 0;
}

port_impl::~port_impl() {

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] X;

	delete[] SIMP;
	delete[] ISIMP;
}

void port_impl::reset() {

	previous_itr_count += itr_tmp;

	itr_tmp = 0;

	// TODO Reset arrays to zero or allocate?
}

void port_impl::add_cols(int n) {

	ASSERT2(A==0,"Already initialized!");

	// FIXME Assumes n==number of constraints, well defined system of equations
	M = n;
	N = 2*M; // at most, equality constraints won't require slack variables

	A = new double[N*M]();
	B = new double[M]();
	C = new double[N]();
	X = new double[N]();

	IS = 2*N;
	SIMP  = new double[IS]();
	ISIMP = new int[IS]();
}

// index[1] ... index[length]
void port_impl::add_eq_row(const int index[], const double value[], int length, double lb, double ub) {

	ASSERT2(M*N>0,"not initialized?");
	ASSERT2(slacks_added<=rows_added&&rows_added<M,"cols_added, rows_added: "<<slacks_added<<", "<<rows_added);

	const int i = rows_added;

	for (int k=1; k<=length; ++k) {

		int j    = index[k] - 1;

		A[j*M+i] = value[k];
	}

	if (lb < ub) {

		int j = M+slacks_added;

		A[j*M+i] = -1;
		B[i] = 0;
		// FIXME Store slack bound in ISIMP and SIMP; make them std::vector
		++slacks_added;
	}
	else {

		B[i] = lb; // == ub
	}

	++rows_added;
}

void port_impl::set_col_bounds(int index, double lb, double ub) {

}

void port_impl::run_simplex() {

}

void port_impl::tighten_col_lb(int i, double& lb) {

}

void port_impl::tighten_col_ub(int i, double& ub) {

}

int port_impl::num_cols() const {
	// FIXME Eliminate as only used for saving dual costs and
	// that should be changes too to save dual costs of rows
}

int port_impl::num_rows() const {

	ASSERT(false);
}

col_status port_impl::col_stat(int i) const {

	ASSERT(false);
}

double port_impl::col_val(int i) const {

}

double port_impl::col_lb(int i) const {

}

double port_impl::col_ub(int i) const {

}

double port_impl::col_dual_val(int i) const {

	ASSERT(false);
}

bool port_impl::is_fixed(int index) const {

	ASSERT(false);
}

void port_impl::dump(const char* file) const {

	ASSERT(false);
}

void port_impl::show_iteration_count() const {

	uint64_t count = previous_itr_count + itr_tmp;

	std::cout << "Simplex iterations: " << count << std::endl;
}

void port_impl::init() {
	// Unused?
}

double port_impl::solve_for(int index, int direction) {

}

}

void cprint_(
		double* A,
		int* M,
		int* N,
		int* IA,
		double* B,
		double* C,
		double* X,
		double* CTX,
		int* IS,
		double* SIMP,
		int* ISIMP,
		int* IE,
		int* ITER,
		int* IPTG,
		int* IAG,
		int* IAS,
		double* U,
		int* DONE)
{
	::itr_tmp = *ITER;
}
