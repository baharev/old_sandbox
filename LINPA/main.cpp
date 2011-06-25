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

extern "C" {

void dlinpr_(double* A, int* M, int* N, int* IA, double* b, double* c, double* x, double* CTX, int* , double* SIMP, int* ISIMP, int* IE);

}

int main() {

	int N = 8;
	int M = 4;

	int IA = 4;
	int IE = 4;

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

	int ISIMP[2*N];
	double SIMP[2*N];

	for (int i=0; i<2*N; ++i) {

		ISIMP[i] = i+1;
		 SIMP[i] = -1;

		ISIMP[i+N] = -(i+1);
		 SIMP[i+N] = 1;
	}

	double c[] = { 12, 2, 10, 14, 0, 1, 2, 0 };

	double x[] = { 0, 0, 0, 0, 0, 0, 0, 0 };

	double CTX;

	int IS = 2*N;

	dlinpr_(A, &M, &N, &IA, b, c, x, &CTX, &IS, SIMP, ISIMP, &IE);

	for (int i=0; i<N; ++i) {
		cout << x[i] << endl;
	}

	return 0;
}
