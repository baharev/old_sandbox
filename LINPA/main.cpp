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

using namespace std;

namespace {
	int itr_tmp = 0;
	int sum_itr = 0;

enum {
	NO_FEASIBLE_SOLUTION = 8
};

}

extern "C" {

void dlinpr_(double* A, int* M, int* N, int* IA, double* b, double* c, double* x, int* MAXITR, double* CTX, int* , double* SIMP, int* ISIMP, int* IE);

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

	itr_tmp = *ITER;

	const int m = *M;
	const int n = *N;

	cout << "===========================================================" << endl;
	cout << "Iteration: " << *ITER << endl;
	cout << "Objective: " << *CTX << endl;
	cout << "X: ";
	for (int i=0; i<n; ++i) {
		cout << X[i] << '\t';
	}
	cout << endl;

	cout << "Residuals of constraints: " << endl;
	for (int i=0; i<m; ++i) {
		double Ax = 0;
		for (int j=0; j<n; ++j) {
			Ax += A[j*m+i]*X[j];
		}
		cout << (i+1) << ":  " << Ax-B[i] << endl;
	}
}

void linpas_(double* A, int* M, int* N, int* IA, double* B, double* C, double* X, int* MAXITR, double* CTX, int* , double* SIMP, int* ISIMP, int* IE, int* ERRCOD);

}

int main() {

	int M = 4;
	int N = 2*M; // If all constraints are double-bounded

	int IA = M;
	int IE = M;

	double A[] = {
			2, 1, 4, 3,
			1, 4, 2, 2,
			3, 2, 3, 4,
			5, 3, 5, 1,
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
	};

	double b[] = { 1, 4, 4, 8 };

	int IS = 2*N;

	int ISIMP[IS];
	double SIMP[IS];

	for (int i=0; i<N; ++i) {

		ISIMP[i] = i+1;
		 SIMP[i] = -1;

		ISIMP[i+N] = -(i+1);
		 SIMP[i+N] = 1;
	}

	double c[] = { 12, 2, 10, 14, 0, 1, 2, 0 };

	//double x[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	//double x[] = { 1, 1, 1, 1, -1, 1, 1, -1 };
	double x[] = { 1, 0.529412, 1, -0.705882, -1, 1, -0.529412, 0.647059 };

	double CTX;

	int MAXITR = 3*N;

	int ERRCOD;

	linpas_(A, &M, &N, &IA, b, c, x, &MAXITR, &CTX, &IS, SIMP, ISIMP, &IE, &ERRCOD);

	sum_itr += itr_tmp;

	itr_tmp = 0;

	cout << "Error code: " << ERRCOD << endl;

	return 0;
}
